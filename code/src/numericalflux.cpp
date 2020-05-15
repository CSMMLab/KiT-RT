#include "numericalflux.h"
#include "upwindflux.h"

NumericalFlux::NumericalFlux( Config* settings ) {}

NumericalFlux* NumericalFlux::Create( Config* settings ) { return new UpwindFlux( settings ); }
