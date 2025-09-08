/* LDCU variant of kEpsilon: applies explicit convection and antidiffusion
   using solver-provided phiLDCUcorr and UstarF. */
#ifndef kEpsilonLdcu_H
#define kEpsilonLdcu_H

#include "kEpsilon.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
class kEpsilonLdcu
:
    public kEpsilon<BasicTurbulenceModel>
{
    using Base = kEpsilon<BasicTurbulenceModel>;

public:
    TypeName("kEpsilonLdcu");

    using typename Base::alphaField;
    using typename Base::rhoField;
    using typename Base::transportModel;

    kEpsilonLdcu
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

