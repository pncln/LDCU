#include "ldcuTurbulenceLib.H"
#include "OSspecific.H"

// A simple function we can reference from a solver so the
// linker includes this library as a dependency and the dynamic
// loader will run its static initializers (registration code).
extern "C" void ldcuTurbulenceEnsureLoaded()
{
    // Intentionally empty. Calling ensures the lib is referenced.
    // Optionally emit a one-time message to aid debugging.
    // Foam::Info<< "libldcuTurbulence: ensureLoaded()" << Foam::nl;
}

