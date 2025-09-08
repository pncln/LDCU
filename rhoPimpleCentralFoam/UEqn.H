// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(U)
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    // Use explicit convection via central-upwind flux
    // Add low-dissipation correction using star velocities
    solve
    (
        UEqn == - fvc::grad(p)
                - fvc::div(phiLDCU, U)
                - fvc::div(phiCorr*(UstarF - fvc::interpolate(U)))
                - fvc::div(momCorrExtra)
    );

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
