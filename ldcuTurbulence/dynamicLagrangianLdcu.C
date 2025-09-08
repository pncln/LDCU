#include "dynamicLagrangianLdcu.H"
#include "fvOptions.H"
#include "upwind.H"

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
dynamicLagrangianLdcu<BasicTurbulenceModel>::dynamicLagrangianLdcu
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
    dynamicLagrangian<BasicTurbulenceModel>
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
void dynamicLagrangianLdcu<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_) return;

    const auto& alpha = this->alpha_;
    const auto& rho   = this->rho_;
    const auto& alphaRhoPhi = this->alphaRhoPhi_;
    const auto& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField S(devSymm(gradU));
    volScalarField magS(mag(S));

    volVectorField Uf(this->filter_(U));
    volSymmTensorField Sf(devSymm(fvc::grad(Uf)));
    volScalarField magSf(mag(Sf));

    volSymmTensorField L(dev(this->filter_(sqr(U)) - (sqr(this->filter_(U)))));
    volSymmTensorField M(2.0*sqr(this->delta())*(this->filter_(magS*S) - 4.0*magSf*Sf));

    volScalarField invT(alpha*rho*(1.0/(this->theta_.value()*this->delta()))*pow(this->flm_*this->fmm_, 1.0/8.0));
    volScalarField LM(L && M);

    // LDCU fields
    const surfaceScalarField& phiCorr =
        this->mesh_.template lookupObject<surfaceScalarField>("phiLDCUcorr");
    const surfaceVectorField& UstarF =
        this->mesh_.template lookupObject<surfaceVectorField>("UstarF");
    const surfaceScalarField phiStar("phiUstar", (UstarF & this->mesh_.Sf()));

    const surfaceScalarField flmFace = fvc::interpolate(this->flm_);
    const surfaceScalarField flmUpw  = upwind<scalar>(this->mesh_, phiStar).interpolate(this->flm_);
    const surfaceScalarField MM(M && M);

    fvScalarMatrix flmEqn
    (
        fvm::ddt(alpha, rho, this->flm_)
      - fvc::div(alphaRhoPhi, this->flm_)
      - fvc::div(phiCorr*(flmUpw - flmFace))
     ==
        invT*LM
      - fvm::Sp(invT, this->flm_)
      + fvOptions(alpha, rho, this->flm_)
    );

    flmEqn.relax();
    fvOptions.constrain(flmEqn);
    flmEqn.solve();
    fvOptions.correct(this->flm_);
    bound(this->flm_, this->flm0_);

    const surfaceScalarField fmmFace = fvc::interpolate(this->fmm_);
    const surfaceScalarField fmmUpw  = upwind<scalar>(this->mesh_, phiStar).interpolate(this->fmm_);

    fvScalarMatrix fmmEqn
    (
        fvm::ddt(alpha, rho, this->fmm_)
      - fvc::div(alphaRhoPhi, this->fmm_)
      - fvc::div(phiCorr*(fmmUpw - fmmFace))
     ==
        invT*MM
      - fvm::Sp(invT, this->fmm_)
      + fvOptions(alpha, rho, this->fmm_)
    );

    fmmEqn.relax();
    fvOptions.constrain(fmmEqn);
    fmmEqn.solve();
    fvOptions.correct(this->fmm_);
    bound(this->fmm_, this->fmm0_);

    this->correctNut(gradU);
}

} // namespace LESModels
} // namespace Foam
