{
    volScalarField& he = thermo.he();

    // Velocity flux for pressure work when he == e
    const surfaceScalarField phivLDCU
    (
        fvc::absolute(phiLDCU/fvc::interpolate(rho), U)
    );

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he)
      + fvc::ddt(rho, K)
      - fvc::div(phiLDCU, he)
      - fvc::div(phiLDCU, K)
      + (
            he.name() == "e"
          ? fvc::div(phivLDCU, p, "div(phiv,p)")
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        fvOptions(rho, he)
    );

    // Low-dissipation kinetic-energy corrections
    const surfaceScalarField KstarF(0.5*magSqr(UstarF));
    // 1) Mass-correction should transport K with U* rather than Uface
    EEqn -= fvc::div(phiCorr*(KstarF - fvc::interpolate(K)));
    // 2) Extra tangential momentum correction contributes mechanical power
    EEqn -= fvc::div(eCorr);

    if (MRF.active())
    {
        EEqn += fvc::div(MRF.phi(), p);
    }

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
