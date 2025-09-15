// LDCU-based continuity equation: use phiLDCU for convective mass flux
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phiLDCU)
     ==
        fvOptions(rho)
    );

    fvOptions.constrain(rhoEqn);
    rhoEqn.solve();
    fvOptions.correct(rho);
}

