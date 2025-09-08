#include "dynamicKEqnLdcu.H"
#include "fvOptions.H"
#include "upwind.H"

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
dynamicKEqnLdcu<BasicTurbulenceModel>::dynamicKEqnLdcu
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
    dynamicKEqn<BasicTurbulenceModel>
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
void dynamicKEqnLdcu<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_) return;

    const auto& alpha = this->alpha_;
    const auto& rho   = this->rho_;
    const auto& alphaRhoPhi = this->alphaRhoPhi_;
    const auto& U = this->U_;
    auto& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volSymmTensorField D(devSymm(tgradU()));
    const volScalarField G(this->GName(), 2.0*nut*(tgradU() && D));
    tgradU.clear();

    volScalarField KK(0.5*(this->filter_(magSqr(U)) - magSqr(this->filter_(U))));
    KK.clamp_min(SMALL);

    // LDCU fields
    const surfaceScalarField& phiCorr =
        this->mesh_.template lookupObject<surfaceScalarField>("phiLDCUcorr");
    const surfaceVectorField& UstarF =
        this->mesh_.template lookupObject<surfaceVectorField>("UstarF");
    const surfaceScalarField phiStar("phiUstar", (UstarF & this->mesh_.Sf()));

    const surfaceScalarField kFace = fvc::interpolate(this->k_);
    const surfaceScalarField kUpw  = upwind<scalar>(this->mesh_, phiStar).interpolate(this->k_);

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, this->k_)
      - fvc::div(alphaRhoPhi, this->k_)
      - fvc::div(phiCorr*(kUpw - kFace))
      - fvm::laplacian(alpha*rho*this->DkEff(), this->k_)
    ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, this->k_)
      - fvm::Sp(this->Ce(D, KK)*alpha*rho*sqrt(this->k_)/this->delta(), this->k_)
      + this->kSource()
      + fvOptions(alpha, rho, this->k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(this->k_);
    bound(this->k_, this->kMin_);

    this->correctNut(D, KK);
}

} // namespace LESModels
} // namespace Foam
