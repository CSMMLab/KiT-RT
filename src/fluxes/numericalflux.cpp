#include "fluxes/numericalflux.hpp"
#include "common/globalconstants.hpp"
#include "fluxes/upwindflux.hpp"

NumericalFluxBase::NumericalFluxBase() {}

NumericalFluxBase* NumericalFluxBase::Create() {
    // TODO: Add Flux options
    return new UpwindFlux();
}
