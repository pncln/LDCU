#include "SpalartAllmarasLdcu.H"
#include "fvOptions.H"
#include "upwind.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
SpalartAllmarasLdcu<BasicTurbulenceModel>::SpalartAllmarasLdcu
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
    SpalartAllmarasBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
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
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


template<class BasicTurbulenceModel>
void SpalartAllmarasLdcu<BasicTurbulenceModel>::correctNut()
{
    // reuse base implementation
    SpalartAllmarasBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut();
}


template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField> SpalartAllmarasLdcu<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    return this->y_;
}


template<class BasicTurbulenceModel>
void SpalartAllmarasLdcu<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_) return;

    const auto& alpha = this->alpha_;
    const auto& rho   = this->rho_;
    const auto& alphaRhoPhi = this->alphaRhoPhi_;
    const auto& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Update model state
    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    const volScalarField ft2(this->ft2(chi));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField dTilda(this->dTilda(chi, fv1, tgradU()));
    volScalarField Stilda(this->Stilda(chi, fv1, tgradU(), dTilda));
    tgradU.clear();

    // LDCU fields
    const surfaceScalarField& phiCorr =
        this->mesh_.template lookupObject<surfaceScalarField>("phiLDCUcorr");
    const surfaceVectorField& UstarF =
        this->mesh_.template lookupObject<surfaceVectorField>("UstarF");
    const surfaceScalarField phiStar("phiUstar", (UstarF & this->mesh_.Sf()));

    const surfaceScalarField nuFace = fvc::interpolate(this->nuTilda_);
    const surfaceScalarField nuUpw  = upwind<scalar>(this->mesh_, phiStar).interpolate(this->nuTilda_);

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, this->nuTilda_)
      - fvc::div(alphaRhoPhi, this->nuTilda_)
      - fvc::div(phiCorr*(nuUpw - nuFace))
      - fvm::laplacian(alpha*rho*this->DnuTildaEff(), this->nuTilda_)
      - this->Cb2_/this->sigmaNut_*alpha()*rho()*magSqr(fvc::grad(this->nuTilda_)()())
     ==
        this->Cb1_*alpha()*rho()*Stilda()*this->nuTilda_()*(scalar(1) - ft2())
      - fvm::Sp
        (
            (this->Cw1_*this->fw(Stilda, dTilda) - this->Cb1_/sqr(this->kappa_)*ft2())
           *alpha()*rho()*this->nuTilda_()/sqr(dTilda()),
            this->nuTilda_
        )
      + fvOptions(alpha, rho, this->nuTilda_)
    );

    nuTildaEqn.ref().relax();
    fvOptions.constrain(nuTildaEqn.ref());
    solve(nuTildaEqn);
    fvOptions.correct(this->nuTilda_);
    bound(this->nuTilda_, dimensionedScalar(this->nuTilda_().dimensions(), Zero));
    this->nuTilda_.correctBoundaryConditions();

    this->correctNut();
}

// runtime selection handled externally; no explicit instantiation

} // namespace RASModels
} // namespace Foam
