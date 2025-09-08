#ifndef realizableKELdcu_H
#define realizableKELdcu_H

#include "realizableKE.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
class realizableKELdcu
:
    public realizableKE<BasicTurbulenceModel>
{
    using Base = realizableKE<BasicTurbulenceModel>;

public:
    TypeName("realizableKELdcu");

    using typename Base::alphaField;
    using typename Base::rhoField;
    using typename Base::transportModel;

    realizableKELdcu
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

