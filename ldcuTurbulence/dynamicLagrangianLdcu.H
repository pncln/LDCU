#ifndef dynamicLagrangianLdcu_H
#define dynamicLagrangianLdcu_H

#include "dynamicLagrangian.H"

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
class dynamicLagrangianLdcu
:
    public dynamicLagrangian<BasicTurbulenceModel>
{
    using Base = dynamicLagrangian<BasicTurbulenceModel>;

public:
    TypeName("dynamicLagrangianLdcu");

    using typename Base::alphaField;
    using typename Base::rhoField;
    using typename Base::transportModel;

    dynamicLagrangianLdcu
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

