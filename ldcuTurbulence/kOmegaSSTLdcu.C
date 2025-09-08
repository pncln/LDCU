#include "kOmegaSSTLdcu.H"
#include "upwind.H"
#include "error.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
kOmegaSSTLdcu<BasicTurbulenceModel>::kOmegaSSTLdcu
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
    kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
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
void kOmegaSSTLdcu<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut(S2);
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaSSTLdcu<BasicTurbulenceModel>::correctNut()
{
    kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaSSTLdcu<BasicTurbulenceModel>::correct()
{
    // replicate base correct(), but swap convection for explicit LDCU
    if (!this->turbulence_)
    {
        return;
    }

    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S2(this->S2(tgradU()));
    volScalarField::Internal GbyNu0(this->GbyNu0(tgradU(), S2));
    volScalarField::Internal G(this->GName(), nut*GbyNu0);

    // Boundary handling per base implementation
    this->omega_.boundaryFieldRef().updateCoeffs();
    this->omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    const volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

    const volScalarField F1(this->F1(CDkOmega));
    const volScalarField F23Field(this->F23());

    // Lookup LDCU fields from registry (guard against missing fields)
    if (!this->mesh_.template foundObject<surfaceScalarField>("phiLDCUcorr")
     || !this->mesh_.template foundObject<surfaceVectorField>("UstarF"))
    {
        FatalErrorInFunction
            << "Required LDCU surface fields 'phiLDCUcorr' and 'UstarF' not found in the mesh registry." << nl
            << "Ensure the solver computes LDCU fluxes before turbulence->correct()." << nl
            << exit(FatalError);
    }

    const surfaceScalarField& phiCorr =
        this->mesh_.template lookupObject<surfaceScalarField>("phiLDCUcorr");
    const surfaceVectorField& UstarF =
        this->mesh_.template lookupObject<surfaceVectorField>("UstarF");
    const surfaceScalarField phiStar("phiUstar", (UstarF & this->mesh_.Sf()));

    // Convenience
    const surfaceScalarField kFace = fvc::interpolate(this->k_);
    const surfaceScalarField omegaFace = fvc::interpolate(this->omega_);
    const surfaceScalarField kUpw = upwind<scalar>(this->mesh_, phiStar).interpolate(this->k_);
    const surfaceScalarField omegaUpw = upwind<scalar>(this->mesh_, phiStar).interpolate(this->omega_);

    {
        const volScalarField::Internal gamma(this->gamma(F1));
        const volScalarField::Internal beta(this->beta(F1));
        GbyNu0 = this->GbyNu(GbyNu0, F23Field, S2());

        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, this->omega_)
          - fvc::div(alphaRhoPhi, this->omega_)
          - fvc::div(phiCorr*(omegaUpw - omegaFace))
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
         ==
            alpha()*rho()*gamma*GbyNu0
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
          - fvm::Sp(alpha()*rho()*beta*this->omega_(), this->omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/this->omega_(),
                this->omega_
            )
          + alpha()*rho()*beta*sqr(this->omegaInf_)
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
          + fvOptions(alpha, rho, this->omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(this->omega_);
        bound(this->omega_, this->omegaMin_);
    }

    {
        tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(alpha, rho, this->k_)
          - fvc::div(alphaRhoPhi, this->k_)
          - fvc::div(phiCorr*(kUpw - kFace))
          - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
         ==
            alpha()*rho()*this->Pk(G)
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
          - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, tgradU()), this->k_)
          + alpha()*rho()*this->betaStar_*this->omegaInf_*this->kInf_
          + this->kSource()
          + fvOptions(alpha, rho, this->k_)
        );

        tgradU.clear();

        kEqn.ref().relax();
        fvOptions.constrain(kEqn.ref());
        solve(kEqn);
        fvOptions.correct(this->k_);
        bound(this->k_, this->kMin_);
    }

    this->correctNut(S2);
}

// Templated class: include instantiations
// No explicit instantiation here; runtime selection handles creation

} // namespace RASModels
} // namespace Foam
