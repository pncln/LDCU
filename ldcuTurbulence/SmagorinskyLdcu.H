#ifndef SmagorinskyLdcu_H
#define SmagorinskyLdcu_H

#include "Smagorinsky.H"

namespace Foam
{
namespace LESModels
{

template<class BasicTurbulenceModel>
class SmagorinskyLdcu
:
    public Smagorinsky<BasicTurbulenceModel>
{
    using Base = Smagorinsky<BasicTurbulenceModel>;
public:
    TypeName("SmagorinskyLdcu");

    using typename Base::alphaField;
    using typename Base::rhoField;
    using typename Base::transportModel;

    SmagorinskyLdcu
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

