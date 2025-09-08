// Ensure LDCU turbulence models are registered even if the dynamic library
// is not loaded via controlDict.libs by compiling the registrations into
// the solver binary.

#include "turbulentFluidThermoModels.H"

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

