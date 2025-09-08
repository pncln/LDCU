#include "realizableKELdcu.H"
#include "fvOptions.H"
#include "upwind.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
realizableKELdcu<BasicTurbulenceModel>::realizableKELdcu
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    realizableKE<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


template<class BasicTurbulenceModel>
void realizableKELdcu<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_) return;

    const auto& alpha = this->alpha_;
    const auto& rho   = this->rho_;
    const auto& alphaRhoPhi = this->alphaRhoPhi_;
    const auto& U = this->U_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    // update nut etc
    realizableKE<BasicTurbulenceModel>::correctNut();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(devSymm(tgradU())));
    volScalarField magS(sqrt(S2));

    volScalarField eta(magS*this->k_/this->epsilon_);
    volScalarField C1(max(eta/(scalar(5) + eta), scalar(0.43)));

    volScalarField G(this->GName(), this->nut_*(tgradU() && devTwoSymm(tgradU())));

    // Boundary updates
    this->epsilon_.boundaryFieldRef().updateCoeffs();
    this->epsilon_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    const surfaceScalarField& phiCorr =
        this->mesh_.template lookupObject<surfaceScalarField>("phiLDCUcorr");
    const surfaceVectorField& UstarF =
        this->mesh_.template lookupObject<surfaceVectorField>("UstarF");
    const surfaceScalarField phiStar("phiUstar", (UstarF & this->mesh_.Sf()));

    const surfaceScalarField epsFace = fvc::interpolate(this->epsilon_);
    const surfaceScalarField epsUpw  = upwind<scalar>(this->mesh_, phiStar).interpolate(this->epsilon_);
    const surfaceScalarField kFace   = fvc::interpolate(this->k_);
    const surfaceScalarField kUpw    = upwind<scalar>(this->mesh_, phiStar).interpolate(this->k_);

    // eps eqn
    // SAF: limit nu to avoid negative viscosity
    volScalarField nuLimited(max(this->nu(), dimensionedScalar(this->nu()().dimensions(), Zero)));

    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, this->epsilon_)
      - fvc::div(alphaRhoPhi, this->epsilon_)
      - fvc::div(phiCorr*(epsUpw - epsFace))
      - fvm::laplacian(alpha*rho*this->DepsilonEff(), this->epsilon_)
     ==
        C1*alpha*rho*magS*this->epsilon_
      - fvm::Sp
        (
            this->C2_*alpha*rho*this->epsilon_/(this->k_ + sqrt(nuLimited*this->epsilon_)),
            this->epsilon_
        )
      + this->epsilonSource()
      + fvOptions(alpha, rho, this->epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(this->epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(this->epsilon_);
    bound(this->epsilon_, this->epsilonMin_);

    // k eqn
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, this->k_)
      - fvc::div(alphaRhoPhi, this->k_)
      - fvc::div(phiCorr*(kUpw - kFace))
      - fvm::laplacian(alpha*rho*this->DkEff(), this->k_)
     ==
        alpha*rho*G
      - fvm::SuSp(2.0/3.0*alpha*rho*divU, this->k_)
      - fvm::Sp(alpha*rho*this->epsilon_/this->k_, this->k_)
      + this->kSource()
      + fvOptions(alpha, rho, this->k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(this->k_);
    bound(this->k_, this->kMin_);

    this->correctNut(tgradU(), S2, magS);
}

} // namespace RASModels
} // namespace Foam
