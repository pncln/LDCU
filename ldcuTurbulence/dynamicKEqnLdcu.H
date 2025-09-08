#ifndef dynamicKEqnLdcu_H
#define dynamicKEqnLdcu_H

#include "dynamicKEqn.H"

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
class dynamicKEqnLdcu
:
    public dynamicKEqn<BasicTurbulenceModel>
{
    using Base = dynamicKEqn<BasicTurbulenceModel>;

public:
    TypeName("dynamicKEqnLdcu");

    using typename Base::alphaField;
    using typename Base::rhoField;
    using typename Base::transportModel;

    dynamicKEqnLdcu
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

