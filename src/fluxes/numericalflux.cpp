#include "fluxes/numericalflux.h"
#include "common/globalconstants.h"
#include "fluxes/upwindflux.h"

NumericalFlux::NumericalFlux() {}

NumericalFlux* NumericalFlux::Create() {
    // TODO: Add Flux options
    return new UpwindFlux();
}
