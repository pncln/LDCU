#ifndef SpalartAllmarasLdcu_H
#define SpalartAllmarasLdcu_H

#include "SpalartAllmaras.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
class SpalartAllmarasLdcu
:
    public SpalartAllmarasBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
{
    using Base = SpalartAllmarasBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>;

protected:
    virtual void correctNut();

    // Provide length scale required by base (same as standard SA)
    virtual tmp<volScalarField> dTilda
    (
        const volScalarField& chi,
        const volScalarField& fv1,
        const volTensorField& gradU
    ) const override;

public:
    using typename Base::alphaField;
    using typename Base::rhoField;
    using typename Base::transportModel;

    TypeName("SpalartAllmarasLdcu");

    SpalartAllmarasLdcu
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    );

    virtual void correct() override;
};

} // namespace RASModels
} // namespace Foam

#endif
