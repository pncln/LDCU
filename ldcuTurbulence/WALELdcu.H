#ifndef WALELdcu_H
#define WALELdcu_H

#include "WALE.H"

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
class WALELdcu
:
    public WALE<BasicTurbulenceModel>
{
    using Base = WALE<BasicTurbulenceModel>;
public:
    TypeName("WALELdcu");

    using typename Base::alphaField;
    using typename Base::rhoField;
    using typename Base::transportModel;

    WALELdcu
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

