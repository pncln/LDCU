#ifndef kEqnLesLdcu_H
#define kEqnLesLdcu_H

#include "kEqn.H"

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
class kEqnLesLdcu
:
    public kEqn<BasicTurbulenceModel>
{
    using Base = kEqn<BasicTurbulenceModel>;

public:
    TypeName("kEqnLesLdcu");

    using typename Base::alphaField;
    using typename Base::rhoField;
    using typename Base::transportModel;

    kEqnLesLdcu
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

} // namespace LESModels
} // namespace Foam

#endif

