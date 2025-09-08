/* Register LDCU turbulence models (RAS/LES) for compressible flows */
#include "turbulentFluidThermoModels.H"
#include "OSspecific.H"

// Emit a message when this library is loaded to confirm registration
namespace {
struct ldcuTurbulenceRegistrationBanner
{
    ldcuTurbulenceRegistrationBanner()
    {
        Foam::Info<< "libldcuTurbulence: registering LDCU RAS/LES models" << Foam::nl;
    }
};

// Static instance constructed at library load time
static const ldcuTurbulenceRegistrationBanner ldcuTurbulenceRegistrationBanner_{};
}

defineTurbulenceModelTypes
(
    geometricOneField,
    volScalarField,
    compressibleTurbulenceModel,
    CompressibleTurbulenceModel,
    ThermalDiffusivity,
    fluidThermo
);

makeBaseTurbulenceModel
(
    geometricOneField,
    volScalarField,
    compressibleTurbulenceModel,
    CompressibleTurbulenceModel,
    ThermalDiffusivity,
    fluidThermo
);

#include "kOmegaSSTLdcu.H"
#include "kEpsilonLdcu.H"
#include "realizableKELdcu.H"
#include "SpalartAllmarasLdcu.H"

#include "kEqnLesLdcu.H"
#include "dynamicKEqnLdcu.H"
#include "dynamicLagrangianLdcu.H"
#include "SmagorinskyLdcu.H"
#include "WALELdcu.H"
#include "sigmaLdcu.H"

// RAS
makeRASModel(kOmegaSSTLdcu);
makeRASModel(kEpsilonLdcu);
makeRASModel(realizableKELdcu);
makeRASModel(SpalartAllmarasLdcu);

// LES
makeLESModel(kEqnLesLdcu);
makeLESModel(dynamicKEqnLdcu);
makeLESModel(dynamicLagrangianLdcu);
makeLESModel(SmagorinskyLdcu);
makeLESModel(WALELdcu);
makeLESModel(sigmaLdcu);
