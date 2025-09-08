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
        F = (sR*FL - sL*FR + sR*sL*(rhoR - rhoL)) / (denom + VSMALL);

        // Low-dissipation correction from LDCU (Part II) — density component
        // Star states (1D-normal) based on local wave speeds
        const scalar rhouL = rhoL*unL;
        const scalar rhouR = rhoR*unR;

    const scalar denomV = (denom + VSMALL);
    const scalar rhoStar = (sR*rhoL - sL*rhoR - (rhouL - rhouR)) / denomV;

    // Momentum star (vector) using full vector flux in n-direction
    const vector FmomL = rhoL*UL*unL + pL*nHat;
    const vector FmomR = rhoR*UR*unR + pR*nHat;
    const vector rhouStarVec = (sR*(rhoL*UL) - sL*(rhoR*UR) - (FmomL - FmomR)) / denomV;
    const scalar rhouStar = rhouStarVec & nHat; // normal component for 1D relations
    const vector Ustar = rhouStarVec / max(rhoStar, VSMALL);
    const scalar uStar = Ustar & nHat;

        const scalar alphaStar = (uStar < 0)
            ? sR/(sR - uStar + VSMALL)
            : sL/(sL - uStar + VSMALL);

    const scalar qRho = minmod
    (
        (uStar - sL)*(rhoStar - rhoR),
        (sR - uStar)*(rhoL - rhoStar)
    );

    const scalar mCorr = alphaStar*qRho; // per unit area [kg/(m^2 s)]
    F += mCorr;

        // Store correction flux and star velocity for later momentum/energy fixes
        const scalar phiCorrInt = mCorr * sArea; // integrated
        phiCorr.primitiveFieldRef()[faceI] = phiCorrInt;
        UstarInt[faceI] = Ustar;

        // Build orthonormal tangential basis (t1,t2)
        vector t1 = vector(0, -nHat.z(), nHat.y());
        if (mag(t1) < SMALL)
        {
            t1 = vector(-nHat.z(), 0, nHat.x());
        }
        if (mag(t1) > VSMALL) t1 /= mag(t1);
        vector t2 = nHat ^ t1; // guaranteed orthogonal

        // Project velocities on tangents
        const scalar vL1 = (UL & t1); const scalar vR1 = (UR & t1);
        const scalar vL2 = (UL & t2); const scalar vR2 = (UR & t2);

        // Star states for rho*v (no pressure in flux), 1: t1, 2: t2
        const scalar rhovL1 = rhoL*vL1; const scalar rhovR1 = rhoR*vR1;
        const scalar rhovL2 = rhoL*vL2; const scalar rhovR2 = rhoR*vR2;

        const scalar rhovStar1 = (sR*rhovL1 - sL*rhovR1 - (rhovL1*unL - rhovR1*unR)) / denom;
        const scalar rhovStar2 = (sR*rhovL2 - sL*rhovR2 - (rhovL2*unL - rhovR2*unR)) / denom;

        const scalar vStar1 = rhovStar1 / max(rhoStar, VSMALL);
        const scalar vStar2 = rhovStar2 / max(rhoStar, VSMALL);

        // a*^+ = sR - u*, a*^- = sL - u*
        const scalar apStar = sR - uStar;
        const scalar amStar = sL - uStar;

        const scalar qRhov1 = minmod((uStar - sL)*(rhovStar1 - rhovR1), (sR - uStar)*(rhovL1 - rhovStar1));
        const scalar qRhov2 = minmod((uStar - sL)*(rhovStar2 - rhovR2), (sR - uStar)*(rhovL2 - rhovStar2));

        // Desired tangential correction flux (integrated)
        vector tangCorr = (alphaStar*qRhov1*sArea)*t1 + (alphaStar*qRhov2*sArea)*t2;

        // Part already accounted by phiCorr * U*_tangential
        const vector UstarTang = Ustar - uStar*nHat;
        vector accounted = phiCorrInt * UstarTang;

        const vector momExtra = tangCorr - accounted;
        momCorrExtra.primitiveFieldRef()[faceI] = momExtra;

        // Exact qE from Part II (2D/3D tangential contributions)
        const scalar rhoPlus1 = max(rhoStar + qRho/(apStar + VSMALL), VSMALL);
        const scalar rhoMinus1 = max(rhoStar + qRho/(amStar + VSMALL), VSMALL);
        const scalar rhovPlus1 = rhovStar1 + qRhov1/(apStar + VSMALL);
        const scalar rhovMinus1 = rhovStar1 + qRhov1/(amStar + VSMALL);

        const scalar rhoPlus2 = max(rhoStar + qRho/(apStar + VSMALL), VSMALL);
        const scalar rhoMinus2 = max(rhoStar + qRho/(amStar + VSMALL), VSMALL);
        const scalar rhovPlus2 = rhovStar2 + qRhov2/(apStar + VSMALL);
        const scalar rhovMinus2 = rhovStar2 + qRhov2/(amStar + VSMALL);

        const scalar tangEnergyDiff = 0.5*
        (
            sqr(rhovPlus1)/rhoPlus1 - sqr(rhovMinus1)/rhoMinus1
          + sqr(rhovPlus2)/rhoPlus2 - sqr(rhovMinus2)/rhoMinus2
        );

        const scalar qE = (apStar*amStar)/(denom) * tangEnergyDiff
                         + 0.5*sqr(uStar)*qRho;

        eCorr.primitiveFieldRef()[faceI] = alphaStar*qE*sArea;
    }

    // Convert to integrated flux (kg/s) with final safety
    const scalar Fint = F * sArea;
    phiInt[faceI] = finite(Fint) ? Fint : phi.primitiveField()[faceI];
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
