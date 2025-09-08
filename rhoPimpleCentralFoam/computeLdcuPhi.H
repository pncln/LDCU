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

// Thermo gamma (ratio of specific heats)
const volScalarField gammaField(thermo.gamma());
const scalarField& gammac = gammaField.primitiveField();

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
label nFallback = 0, nNearChar = 0, nClamped = 0, nNonFinite = 0;
// Start-up ramp for corrections to avoid shocks at t=0 (first ~200 time steps)
const scalar ldcuRamp = min(scalar(1), scalar(runTime.timeIndex())/scalar(200));
scalar minRhoStarDiag = GREAT, minPStarDiag = GREAT;
vectorField& UstarInt = UstarF.primitiveFieldRef();
for (label faceI = 0; faceI < nIntFaces; ++faceI)
{
    const label cl = owner[faceI];
    const label cr = neigh[faceI];

    const vector& Sf_i = Sf[faceI];
    const scalar sArea = magSf[faceI] + VSMALL;
    const vector nHat = Sf_i / sArea;

    const vector& UL = Uc[cl];
    const vector& UR = Uc[cr];
    const scalar rhoL = max(rhoc[cl], VSMALL);
    const scalar rhoR = max(rhoc[cr], VSMALL);
    const scalar pL = max(pc[cl], VSMALL);
    const scalar pR = max(pc[cr], VSMALL);
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
        const scalar uStar = Ustar & nHat;

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
        const scalar EL = pL/(gammaFace - 1.0) + 0.5*rhoL*(UL & UL);
        const scalar ER = pR/(gammaFace - 1.0) + 0.5*rhoR*(UR & UR);
        const scalar f4L = unL*(EL + pL);
        const scalar f4R = unR*(ER + pR);
        const scalar EStar = (sR*EL - sL*ER - (f4R - f4L)) / denomV;
        const scalar pStar = (gammaFace - 1.0) * (EStar - 0.5*magSqr(rhouStarVec)/max(rhoStar, VSMALL));
        minRhoStarDiag = min(minRhoStarDiag, rhoStar);
        minPStarDiag = min(minPStarDiag, pStar);
        const scalar dp = max(pStar - 1e-13, 0.0);

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
        dE += 0.5*
            ((sqr(mStar1 - dm1/(ap1 + VSMALL)) / max(rhoStar - drho/(ap1 + VSMALL), VSMALL)
             - sqr(mStar1 - dm1/(am1 + VSMALL)) / max(rhoStar - drho/(am1 + VSMALL), VSMALL))
            ) / (denom + VSMALL);
        dE += 0.5*
            ((sqr(mStar2 - dm2/(ap1 + VSMALL)) / max(rhoStar - drho/(ap1 + VSMALL), VSMALL)
             - sqr(mStar2 - dm2/(am1 + VSMALL)) / max(rhoStar - drho/(am1 + VSMALL), VSMALL))
            ) / (denom + VSMALL);

        // Scale as in reference (guard against near-characteristic speeds)
        const bool pos = (uStar > 0);
        const scalar waveTol = 1e-8*max(mag(sR), mag(sL)) + VSMALL;
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
            const scalar scaleCap = 2.0;
            scaleR = Foam::max(Foam::min(scaleR, scaleCap), -scaleCap);
            drho *= scaleR;
            dm1  *= scaleR;
            dm2  *= scaleR;
            dE   *= scaleE;
        }

        // Zero if u* exceeds wave speeds
        if (mag(uStar) >= max(sR, -sL)) { drho = 0; dm1 = 0; dm2 = 0; dE = 0; }

        // Cap mass correction vs baseline flux magnitude
        const scalar cap = 0.25*(mag(Fbase) + SMALL);
        if (mag(drho) > cap)
        {
            drho = copysign(cap, drho);
            ++nClamped;
        }

        // Add kinetic contribution
        dE += 0.5*sqr(uStar)*drho;
        // Apply start-up ramp to all corrections
        drho *= ldcuRamp;
        dm1  *= ldcuRamp;
        dm2  *= ldcuRamp;
        dE   *= ldcuRamp;

        // Ensure finite values; otherwise zero the corrections
        if (!finite(drho)) { drho = 0.0; ++nNonFinite; }
        if (!finite(dm1)) dm1 = 0.0;
        if (!finite(dm2)) dm2 = 0.0;
        if (!finite(dE))  { dE  = 0.0; ++nNonFinite; }
        if (!finite(uStar)) { Ustar = (UL+UR)/2.0; }

        // Corrected mass flux per unit area
        F -= drho;

        // Store integrated corrections for momentum/energy equations
        const scalar phiCorrInt = drho * sArea; // integrated mass correction
        phiCorr.primitiveFieldRef()[faceI] = phiCorrInt;
        UstarInt[faceI] = Ustar;

        // Tangential momentum corrections: split accounted vs extra
        const vector tangCorr = (dm1*sArea)*t1 + (dm2*sArea)*t2;
        const vector UstarTang = Ustar - uStar*nHat;
        const vector accounted = phiCorrInt * UstarTang;
        const vector momExtra = tangCorr - accounted;
        momCorrExtra.primitiveFieldRef()[faceI] =
            vector
            (
                finite(momExtra.x()) ? momExtra.x() : 0.0,
                finite(momExtra.y()) ? momExtra.y() : 0.0,
                finite(momExtra.z()) ? momExtra.z() : 0.0
            );

        // Energy correction integrated
        {
            const scalar eInt = dE*sArea;
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
