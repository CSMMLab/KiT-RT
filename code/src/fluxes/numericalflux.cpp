#include "fluxes/numericalflux.h"
#include "fluxes/upwindflux.h"
#include "settings/globalconstants.h"

NumericalFlux::NumericalFlux() {}

NumericalFlux* NumericalFlux::Create() {
    // TODO: Add Flux options
    return new UpwindFlux();
}
