#ifndef sigmaLdcu_H
#define sigmaLdcu_H

#include "sigma.H"

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
class sigmaLdcu
:
    public sigma<BasicTurbulenceModel>
{
    using Base = sigma<BasicTurbulenceModel>;
public:
    TypeName("sigmaLdcu");

    using typename Base::alphaField;
    using typename Base::rhoField;
    using typename Base::transportModel;

    sigmaLdcu
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    ) : Base(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName, type) {}
};

} // namespace LESModels
} // namespace Foam

#endif

