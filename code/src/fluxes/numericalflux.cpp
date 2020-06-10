#include "fluxes/numericalflux.h"
#include "fluxes/upwindflux.h"

NumericalFlux::NumericalFlux( Config* settings ) {}

NumericalFlux* NumericalFlux::Create( Config* settings ) {
    // TODO: Add Flux options
    return new UpwindFlux( settings );
}
