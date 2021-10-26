#include "fluxes/numericalflux.hpp"
#include "common/globalconstants.hpp"
#include "fluxes/upwindflux.hpp"

NumericalFlux::NumericalFlux() {}

NumericalFlux* NumericalFlux::Create() {
    // TODO: Add Flux options
    return new UpwindFlux();
}
