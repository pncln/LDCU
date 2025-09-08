/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    LDCU variant of kOmegaSST that applies explicit LDCU-style convection
    and antidiffusion to k and omega using phiLDCU and phiLDCUcorr/UstarF
    surface fields produced by the solver.
-------------------------------------------------------------------------------*/

#ifndef kOmegaSSTLdcu_H
#define kOmegaSSTLdcu_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "kOmegaSSTBase.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
class kOmegaSSTLdcu
:
    public kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
{
    // Private Member Functions
        kOmegaSSTLdcu(const kOmegaSSTLdcu&) = delete;
        void operator=(const kOmegaSSTLdcu&) = delete;

protected:
        virtual void correctNut(const volScalarField& S2);
        virtual void correctNut();

public:
    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("kOmegaSSTLdcu");

    // Constructors
    kOmegaSSTLdcu
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

    //- Override model correct to apply LDCU-style convection/antidiffusion
    virtual void correct();
};

} // End namespace RASModels
} // End namespace Foam

// (NoRepository instantiation not used here; compiled via library object)

#endif

// ************************************************************************* //
