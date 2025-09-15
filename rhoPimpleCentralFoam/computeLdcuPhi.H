// Compute central-upwind (Kurganov-type) mass flux phiLDCU for use
// in convective terms. This overrides fvSchemes convection implicitly
// by providing an explicit mass flux based on left/right states and
// local wave speeds. Boundary fluxes fall back to existing phi.

#include <cmath>

// Create phiLDCU by copying topology/boundaries from phi
surfaceScalarField phiLDCU("phiLDCU", phi);
// Mass-correction part (alpha* q^rho), integrated [kg/s]
surfaceScalarField phiCorr
(
    IOobject
    (
        "phiLDCUcorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", phi.dimensions(), 0.0)
);
// Star velocity at faces used for low-dissipation correction
surfaceVectorField UstarF
(
    IOobject
    (
        "UstarF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(U)
);

// Additional tangential momentum correction beyond q^rho * U*
surfaceVectorField momCorrExtra
(
    IOobject
    (
        "momCorrExtra",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", phi.dimensions()*dimVelocity, vector::zero)
);

// Energy correction associated with extra tangential momentum transport
surfaceScalarField eCorr
(
    IOobject
    (
        "eCorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", phi.dimensions()*dimVelocity*dimVelocity, 0.0)
);

// Access mesh data
const labelUList& owner = mesh.owner();
const labelUList& neigh = mesh.neighbour();
const vectorField& Sf = mesh.Sf();
const scalarField magSf(mag(Sf));

// Internal fields
const vectorField& Uc = U.primitiveField();
const scalarField& rhoc = rho.primitiveField();
const scalarField& pc = p.primitiveField();

// Cell and face centres
const pointField& Cc = mesh.C();
const pointField& Cf = mesh.Cf();

// Thermo gamma (ratio of specific heats)
const volScalarField gammaField(thermo.gamma());
const scalarField& gammac = gammaField.primitiveField();

// Conservative variables for reconstruction
tmp<volVectorField> tRhoU = rho*U;
const volVectorField& rhoUField = tRhoU();
tmp<volScalarField> tEtot = p/(gammaField - scalar(1)) + 0.5*magSqr(rhoUField)/rho;
const volScalarField& Etot = tEtot();

// Gradients for MUSCL-type reconstruction on conservative vars
tmp<volVectorField> tGradRho = fvc::grad(rho);
const volVectorField& gradRho = tGradRho();
tmp<volTensorField> tGradRhoU = fvc::grad(rhoUField);
const volTensorField& gradRhoU = tGradRhoU();
tmp<volVectorField> tGradE = fvc::grad(Etot);
const volVectorField& gradE = tGradE();

// Precompute neighbour-based bounds (Barth-Jespersen stencil)
const label nCells = mesh.nCells();
const vectorField& mC = rhoUField.primitiveField();
const scalarField& Ec = Etot.primitiveField();
scalarField rhoMin(nCells), rhoMax(nCells), EMin(nCells), EMax(nCells);
scalarField mxMin(nCells), mxMax(nCells), myMin(nCells), myMax(nCells), mzMin(nCells), mzMax(nCells);
forAll(rhoc, cellI)
{
    const vector& mc = mC[cellI];
    rhoMin[cellI] = rhoc[cellI];
    rhoMax[cellI] = rhoc[cellI];
    EMin[cellI] = Ec[cellI];
    EMax[cellI] = Ec[cellI];
    mxMin[cellI] = mc.x(); mxMax[cellI] = mc.x();
    myMin[cellI] = mc.y(); myMax[cellI] = mc.y();
    mzMin[cellI] = mc.z(); mzMax[cellI] = mc.z();
}
for (label faceI = 0; faceI < owner.size(); ++faceI)
{
    const label cl = owner[faceI];
    const label cr = neigh[faceI];
    const scalar rhcr = rhoc[cr];
    const scalar rhcl = rhoc[cl];
    rhoMin[cl] = min(rhoMin[cl], rhcr); rhoMax[cl] = max(rhoMax[cl], rhcr);
    rhoMin[cr] = min(rhoMin[cr], rhcl); rhoMax[cr] = max(rhoMax[cr], rhcl);
    const scalar Ecr = Ec[cr]; const scalar Ecl = Ec[cl];
    EMin[cl] = min(EMin[cl], Ecr); EMax[cl] = max(EMax[cl], Ecr);
    EMin[cr] = min(EMin[cr], Ecl); EMax[cr] = max(EMax[cr], Ecl);
    const vector& mcr = mC[cr]; const vector& mcl = mC[cl];
    mxMin[cl] = min(mxMin[cl], mcr.x()); mxMax[cl] = max(mxMax[cl], mcr.x());
    myMin[cl] = min(myMin[cl], mcr.y()); myMax[cl] = max(myMax[cl], mcr.y());
    mzMin[cl] = min(mzMin[cl], mcr.z()); mzMax[cl] = max(mzMax[cl], mcr.z());
    mxMin[cr] = min(mxMin[cr], mcl.x()); mxMax[cr] = max(mxMax[cr], mcl.x());
    myMin[cr] = min(myMin[cr], mcl.y()); myMax[cr] = max(myMax[cr], mcl.y());
    mzMin[cr] = min(mzMin[cr], mcl.z()); mzMax[cr] = max(mzMax[cr], mcl.z());
}

// Venkatakrishnan smooth limiter for scalars
auto venkat = [](const scalar qc, const scalar qmin, const scalar qmax, const scalar dq, const scalar eps2) -> scalar
{
    if (dq > 0)
    {
        const scalar a = qmax - qc;
        const scalar num = a*a + 2.0*a*dq + eps2;
        const scalar den = a*a + a*dq + 2.0*dq*dq + eps2 + VSMALL;
        return max(scalar(0), min(scalar(1), num/den));
    }
    else if (dq < 0)
    {
        const scalar b = qmin - qc; // negative or zero
        const scalar num = b*b + 2.0*b*dq + eps2;
        const scalar den = b*b + b*dq + 2.0*dq*dq + eps2 + VSMALL;
        return max(scalar(0), min(scalar(1), num/den));
    }
    return scalar(1);
};


// Helper: minmod limiter
auto minmod = [](const scalar a, const scalar b) -> scalar
{
    if (a*b <= 0) return 0.0;
    return (mag(a) < mag(b)) ? a : b;
};

auto finite = [](const scalar x) -> bool { return std::isfinite(double(x)); };

// Internal face fluxes
scalarField& phiInt = phiLDCU.primitiveFieldRef();

const label nIntFaces = owner.size();
// Precompute face-interpolated U for safe fallback use
tmp<surfaceVectorField> tUface = fvc::interpolate(U);
const surfaceVectorField& Uface = tUface();
label nFallback = 0, nNearChar = 0, nClamped = 0, nNonFinite = 0;
// Start-up ramp for corrections to avoid shocks at t=0 (first ~200 time steps)
const scalar ldcuRamp = min(scalar(1), scalar(runTime.timeIndex())/scalar(200));
// Tangential corrections (dm, dE) disabled for an initial window,
// then ramped in. Configurable via controlDict:
//   ldcuTangentialOffSteps (default 5000)
//   ldcuTangentialRampSteps (default 3000)
const label tangentialOffSteps =
    runTime.controlDict().lookupOrDefault<label>("ldcuTangentialOffSteps", 5000);
const label tangentialRampSteps =
    runTime.controlDict().lookupOrDefault<label>("ldcuTangentialRampSteps", 3000);
const bool tangentialBaseOn = (runTime.timeIndex() >= tangentialOffSteps);
const scalar tangRamp =
    tangentialBaseOn
  ? min(scalar(1), scalar(runTime.timeIndex() - tangentialOffSteps)/scalar(tangentialRampSteps))
  : scalar(0);
scalar minRhoStarDiag = GREAT, minPStarDiag = GREAT;
vectorField& UstarInt = UstarF.primitiveFieldRef();
for (label faceI = 0; faceI < nIntFaces; ++faceI)
{
    const label cl = owner[faceI];
    const label cr = neigh[faceI];

    const vector& Sf_i = Sf[faceI];
    const scalar sArea = magSf[faceI] + VSMALL;
    const vector nHat = Sf_i / sArea;

    // Reconstruct left/right conservative states at face center (BJ limiter)
    const scalar rhoLc = rhoc[cl];
    const scalar rhoRc = rhoc[cr];
    const vector mLc = mC[cl];
    const vector mRc = mC[cr];
    const scalar ELc = Ec[cl];
    const scalar ERc = Ec[cr];

    const vector dL = Cf[faceI] - Cc[cl];
    const vector dR = Cf[faceI] - Cc[cr];

    // Unlimited reconstructions of conservative variables
    const scalar dRhoL = (gradRho[cl] & dL);
    const scalar dRhoR = (gradRho[cr] & dR);
    const vector dML = (gradRhoU[cl] & dL);
    const vector dMR = (gradRhoU[cr] & dR);
    const scalar dEL = (gradE[cl] & dL);
    const scalar dER = (gradE[cr] & dR);

    // Venkatakrishnan limiter factors (use local variable ranges for epsilon)
    const scalar kappa = 0.3; // tunable, smoothness parameter
    const scalar eps2_rhoL = sqr(kappa*(rhoMax[cl] - rhoMin[cl]) + VSMALL);
    const scalar eps2_rhoR = sqr(kappa*(rhoMax[cr] - rhoMin[cr]) + VSMALL);
    const scalar eps2_EL = sqr(kappa*(EMax[cl] - EMin[cl]) + VSMALL);
    const scalar eps2_ER = sqr(kappa*(EMax[cr] - EMin[cr]) + VSMALL);
    const scalar eps2_mxL = sqr(kappa*(mxMax[cl] - mxMin[cl]) + VSMALL);
    const scalar eps2_myL = sqr(kappa*(myMax[cl] - myMin[cl]) + VSMALL);
    const scalar eps2_mzL = sqr(kappa*(mzMax[cl] - mzMin[cl]) + VSMALL);
    const scalar eps2_mxR = sqr(kappa*(mxMax[cr] - mxMin[cr]) + VSMALL);
    const scalar eps2_myR = sqr(kappa*(myMax[cr] - myMin[cr]) + VSMALL);
    const scalar eps2_mzR = sqr(kappa*(mzMax[cr] - mzMin[cr]) + VSMALL);

    const scalar phiRhoL = venkat(rhoLc, rhoMin[cl], rhoMax[cl], dRhoL, eps2_rhoL);
    const scalar phiRhoR = venkat(rhoRc, rhoMin[cr], rhoMax[cr], dRhoR, eps2_rhoR);
    const scalar phiEL = venkat(ELc, EMin[cl], EMax[cl], dEL, eps2_EL);
    const scalar phiER = venkat(ERc, EMin[cr], EMax[cr], dER, eps2_ER);
    const scalar phiMxL = venkat(mLc.x(), mxMin[cl], mxMax[cl], dML.x(), eps2_mxL);
    const scalar phiMyL = venkat(mLc.y(), myMin[cl], myMax[cl], dML.y(), eps2_myL);
    const scalar phiMzL = venkat(mLc.z(), mzMin[cl], mzMax[cl], dML.z(), eps2_mzL);
    const scalar phiMxR = venkat(mRc.x(), mxMin[cr], mxMax[cr], dMR.x(), eps2_mxR);
    const scalar phiMyR = venkat(mRc.y(), myMin[cr], myMax[cr], dMR.y(), eps2_myR);
    const scalar phiMzR = venkat(mRc.z(), mzMin[cr], mzMax[cr], dMR.z(), eps2_mzR);

    const scalar rhoL = max(rhoLc + phiRhoL*dRhoL, VSMALL);
    const scalar rhoR = max(rhoRc + phiRhoR*dRhoR, VSMALL);
    const vector mL(mLc.x() + phiMxL*dML.x(), mLc.y() + phiMyL*dML.y(), mLc.z() + phiMzL*dML.z());
    const vector mR(mRc.x() + phiMxR*dMR.x(), mRc.y() + phiMyR*dMR.y(), mRc.z() + phiMzR*dMR.z());
    const scalar EL = max(ELc + phiEL*dEL, VSMALL);
    const scalar ER = max(ERc + phiER*dER, VSMALL);

    // Primitive from conservative reconstructions
    const vector UL = mL/max(rhoL, VSMALL);
    const vector UR = mR/max(rhoR, VSMALL);
    const scalar pL = max((gammac[cl] - 1.0)*(EL - 0.5*magSqr(mL)/max(rhoL, VSMALL)), VSMALL);
    const scalar pR = max((gammac[cr] - 1.0)*(ER - 0.5*magSqr(mR)/max(rhoR, VSMALL)), VSMALL);
    const scalar gammaL = max(gammac[cl], SMALL);
    const scalar gammaR = max(gammac[cr], SMALL);

    const scalar unL = (UL & nHat);
    const scalar unR = (UR & nHat);

    const scalar aL = Foam::sqrt(max(gammaL*pL/rhoL, SMALL));
    const scalar aR = Foam::sqrt(max(gammaR*pR/rhoR, SMALL));

    // Wave speeds (Kurganov), clamped with zero as in reference code
    const scalar sL = min(min(unL - aL, unR - aR), scalar(0));
    const scalar sR = max(max(unL + aL, unR + aR), scalar(0));

    const scalar denom = (sR - sL);
    if (!finite(sL) || !finite(sR) || mag(denom) < SMALL)
    {
        // Degenerate wave speeds or bad numerics: fall back to central average
        const scalar fL = rhoL*unL;
        const scalar fR = rhoR*unR;
        phiInt[faceI] = 0.5*(fL + fR) * sArea;
        // Keep default zeros for corrections and initial interpolated Ustar
        ++nFallback;
        continue;
    }

    // Left/right mass flux per unit area
    const scalar FL = rhoL*unL;
    const scalar FR = rhoR*unR;

    scalar F = 0.0;
    if (sL >= 0)
    {
        F = FL;
    }
    else if (sR <= 0)
    {
        F = FR;
    }
    else
    {
        // Central-upwind baseline flux (Kurganov-Tadmor type)
        const scalar Fbase = (sR*FL - sL*FR + sR*sL*(rhoR - rhoL)) / (denom + VSMALL);
        F = Fbase;

        // Star states (1D-normal) based on local wave speeds
        const scalar rhouL = rhoL*unL;
        const scalar rhouR = rhoR*unR;

    const scalar denomV = (denom + VSMALL);
        const scalar rhoStar = (sR*rhoL - sL*rhoR - (rhouL - rhouR)) / denomV;

        // Momentum star (vector) using full vector flux in n-direction
        const vector FmomL = rhoL*UL*unL + pL*nHat;
        const vector FmomR = rhoR*UR*unR + pR*nHat;
        const vector rhouStarVec = (sR*(rhoL*UL) - sL*(rhoR*UR) - (FmomL - FmomR)) / denomV;
        vector Ustar = rhouStarVec / max(rhoStar, VSMALL);
        scalar uStar = Ustar & nHat;

        // Build orthonormal tangential basis (t1,t2)
        vector t1 = vector(0, -nHat.z(), nHat.y());
        if (mag(t1) < SMALL) { t1 = vector(-nHat.z(), 0, nHat.x()); }
        if (mag(t1) > VSMALL) t1 /= mag(t1);
        vector t2 = nHat ^ t1; // orthogonal

        // Tangential momenta
        const scalar vL1 = (UL & t1), vR1 = (UR & t1);
        const scalar vL2 = (UL & t2), vR2 = (UR & t2);
        const scalar mL1 = rhoL*vL1, mR1 = rhoR*vR1;
        const scalar mL2 = rhoL*vL2, mR2 = rhoR*vR2;
        const scalar mStar1 = rhouStarVec & t1;
        const scalar mStar2 = rhouStarVec & t2;

        // Energy star (perfect gas)
        const scalar gammaFace = 0.5*(gammaL + gammaR);
        // Explicit physical energy from p and U for fluxes
        const scalar EphysL = pL/(gammaFace - 1.0) + 0.5*rhoL*(UL & UL);
        const scalar EphysR = pR/(gammaFace - 1.0) + 0.5*rhoR*(UR & UR);
        const scalar f4L = unL*(EphysL + pL);
        const scalar f4R = unR*(EphysR + pR);
        const scalar EStar = (sR*EphysL - sL*EphysR - (f4R - f4L)) / denomV;
        const scalar pStar = (gammaFace - 1.0) * (EStar - 0.5*magSqr(rhouStarVec)/max(rhoStar, VSMALL));
        minRhoStarDiag = min(minRhoStarDiag, rhoStar);
        minPStarDiag = min(minPStarDiag, pStar);
        const scalar dp = max(pStar - 1e-13, 0.0);
        // Simple shock sensor based on relative pressure jump
        const scalar pAvg = 0.5*(pL + pR) + SMALL;
        const scalar shockRef = 0.1;
        scalar shockS = min(1.0, mag(pR - pL)/(shockRef*pAvg));
        shockS = max(0.0, shockS);
        // Additional sensor for strongly negative pStar (pre-shock rarefaction)
        scalar starShock = 0.0;
        if (pStar < 0)
        {
            starShock = min(1.0, (-pStar)/(shockRef*pAvg));
        }

        // Optional positivity-preserving fallback (local Lax-Friedrichs)
        // If star states are non-physical, switch to LLF flux and suppress corrections
        const scalar aMax = max(mag(unL) + aL, mag(unR) + aR);
        if (pStar <= 0 || rhoStar <= VSMALL)
        {
            const scalar Fllf = 0.5*(FL + FR) - 0.5*aMax*(rhoR - rhoL);
            const scalar FintLLF = Fllf * sArea;
            phiInt[faceI] = finite(FintLLF) ? FintLLF : 0.5*(FL + FR)*sArea;
            // Suppress all LDCU corrections on this face
            phiCorr.primitiveFieldRef()[faceI] = 0.0;
            UstarInt[faceI] = Uface[faceI];
            momCorrExtra.primitiveFieldRef()[faceI] = vector::zero;
            eCorr.primitiveFieldRef()[faceI] = 0.0;
            continue;
        }

        // a*^+ = sR - u*, a*^- = sL - u*
        const scalar ap1 = sR - uStar;
        const scalar am1 = sL - uStar;

        // Density/tangential corrections with limiter
        scalar drho = minmod(-am1*(rhoL - rhoStar), ap1*(rhoStar - rhoR));
        scalar dm1 = minmod(-am1*(mL1 - mStar1), ap1*(mStar1 - mR1));
        scalar dm2 = minmod(-am1*(mL2 - mStar2), ap1*(mStar2 - mR2));

        const scalar rhoLmod = am1*rhoStar - drho;
        const scalar rhoRmod = ap1*rhoStar - drho;
        const scalar vStar1 = Ustar & t1;
        const scalar vStar2 = Ustar & t2;
        const scalar dmmax =
            (dp > VSMALL && rhoStar > VSMALL)
          ? Foam::sqrt(max(0.0, -2.0*rhoLmod*rhoRmod*dp/(rhoStar*(gammaFace-1.0))))
          : 0.0;
        if (dmmax > 0)
        {
            const scalar diff1 = dm1 - vStar1*drho;
            if (mag(diff1) > dmmax) dm1 = vStar1*drho + std::copysign(dmmax, diff1);
            const scalar diff2 = dm2 - vStar2*drho;
            if (mag(diff2) > dmmax) dm2 = vStar2*drho + std::copysign(dmmax, diff2);
        }

        // Energy correction (sum over tangential components)
        scalar dE = 0.0;
        const bool tangentialEnabled = tangentialBaseOn && (pStar >= 0);
        if (tangentialEnabled)
        {
            // additional safety for energy correction denominators
            const scalar rhoMinSafe = 1e-12;
            const scalar denAp = max(rhoStar - drho/(ap1 + VSMALL), rhoMinSafe);
            const scalar denAm = max(rhoStar - drho/(am1 + VSMALL), rhoMinSafe);

            dE += 0.5*
                ((sqr(mStar1 - dm1/(ap1 + VSMALL)) / denAp
                 - sqr(mStar1 - dm1/(am1 + VSMALL)) / denAm)
                ) / (denom + VSMALL);
            dE += 0.5*
                ((sqr(mStar2 - dm2/(ap1 + VSMALL)) / denAp
                 - sqr(mStar2 - dm2/(am1 + VSMALL)) / denAm)
                ) / (denom + VSMALL);
        }

        // Scale as in reference (guard against near-characteristic speeds)
        const bool pos = (uStar > 0);
        const scalar waveTol = 1e-6*(mag(sR) + mag(sL) + mag(aL) + mag(aR) + mag(unL) + mag(unR)) + SMALL;
        if (mag(ap1) < waveTol || mag(am1) < waveTol)
        {
            drho = 0; dm1 = 0; dm2 = 0; dE = 0;
            ++nNearChar;
        }
        else
        {
            scalar scaleR = pos ? (sL/(am1)) : (sR/(ap1));
            const scalar scaleE = pos ? (sL*ap1) : (sR*am1);
            // Clamp scale to avoid excessive amplification
            const scalar scaleCap = 0.3;
            scaleR = Foam::max(Foam::min(scaleR, scaleCap), -scaleCap);
            drho *= scaleR;
            if (tangentialEnabled)
            {
                dm1  *= scaleR;
                dm2  *= scaleR;
                dE   *= scaleE;
            }
            else
            {
                dm1 = 0; dm2 = 0; dE = 0;
            }
        }

        // Zero if u* exceeds wave speeds
        if (mag(uStar) >= max(sR, -sL)) { drho = 0; dm1 = 0; dm2 = 0; dE = 0; }

        // Cap mass correction vs baseline flux magnitude (more conservative early on)
        const scalar clampFactor = 0.02 + 0.10*min(scalar(1), scalar(runTime.timeIndex())/scalar(200));
        const scalar cap = clampFactor*(mag(Fbase) + SMALL);
        if (mag(drho) > cap)
        {
            drho = copysign(cap, drho);
            ++nClamped;
        }

        // If p* is negative, suppress all corrections on this face
        if (pStar <= 0) { drho = 0; dm1 = 0; dm2 = 0; dE = 0; }

        // Add kinetic contribution
        dE += 0.5*sqr(uStar)*drho;
        // Tangential ramp-in after the start window
        dm1 *= tangRamp;
        dm2 *= tangRamp;
        dE  *= tangRamp;
        // Apply start-up ramp and shock blending (quadratic) to all corrections
        const scalar blend = ldcuRamp * sqr(1.0 - shockS) * sqr(1.0 - starShock);
        drho *= blend;
        dm1  *= blend;
        dm2  *= blend;
        dE   *= blend;

        // Ensure finite values; otherwise zero the corrections
        if (!finite(drho)) { drho = 0.0; ++nNonFinite; }
        if (!finite(dm1)) dm1 = 0.0;
        if (!finite(dm2)) dm2 = 0.0;
        if (!finite(dE))  { dE  = 0.0; ++nNonFinite; }
        if (!finite(uStar))
        {
            Ustar = (UL+UR)/2.0;
            uStar = Ustar & nHat; // refresh normal projection
        }

        // Corrected mass flux per unit area
        F -= drho;

        // Store integrated corrections for momentum/energy equations
        scalar phiCorrInt = drho * sArea; // integrated mass correction
        // In strict drho-only mode (tangential disabled), suppress phiCorr usage
        if (!tangentialEnabled)
        {
            phiCorrInt = 0.0;
        }
        phiCorr.primitiveFieldRef()[faceI] = phiCorrInt;
        UstarInt[faceI] = Ustar;

        // Tangential momentum corrections: split accounted vs extra
        // Convert tangential momentum density correction to a flux using the
        // advecting normal speed u*: ensures units [kg*m/s^2].
        const vector tangCorrFlux = uStar * ((dm1*sArea)*t1 + (dm2*sArea)*t2);
        const vector UstarTang   = Ustar - uStar*nHat;
        const vector accounted   = phiCorrInt * UstarTang;
        const vector momExtra    = tangCorrFlux - accounted;
        momCorrExtra.primitiveFieldRef()[faceI] =
            vector
            (
                finite(momExtra.x()) ? momExtra.x() : 0.0,
                finite(momExtra.y()) ? momExtra.y() : 0.0,
                finite(momExtra.z()) ? momExtra.z() : 0.0
            );

        // Energy correction integrated
        {
            // In drho-only mode, suppress energy correction
            const scalar eInt = tangentialEnabled ? dE*sArea : 0.0;
            eCorr.primitiveFieldRef()[faceI] = finite(eInt) ? eInt : 0.0;
        }
    }

    // Convert to integrated flux (kg/s) with final safety
    const scalar Fint = F * sArea;
    if (!finite(Fint)) { phiInt[faceI] = 0.5*(FL + FR)*sArea; ++nNonFinite; }
    else { phiInt[faceI] = Fint; }
}

// For boundaries, fall back to the existing mass flux
forAll(phiLDCU.boundaryField(), patchI)
{
    phiLDCU.boundaryFieldRef()[patchI] = phi.boundaryField()[patchI];
    // No correction and use interpolated velocity on boundaries
    phiCorr.boundaryFieldRef()[patchI] = 0.0;
    UstarF.boundaryFieldRef()[patchI] = fvc::interpolate(U)().boundaryField()[patchI];
    momCorrExtra.boundaryFieldRef()[patchI] = vector::zero;
    eCorr.boundaryFieldRef()[patchI] = 0.0;
}

// Diagnostics (first few timesteps) – print if any events occurred
if (nFallback || nNearChar || nClamped || nNonFinite)
{
    Info<< "LDCU diagnostics: fallback=" << nFallback
        << " nearChar=" << nNearChar
        << " clamped=" << nClamped
        << " nonFinite=" << nNonFinite
        << " min(rhoStar)=" << minRhoStarDiag
        << " min(pStar)=" << minPStarDiag << nl;
}
