#include "numericalflux.h"
#include "upwindflux.h"

NumericalFlux::NumericalFlux( Settings* settings ) {}

NumericalFlux* NumericalFlux::Create( Settings* settings ) { return new UpwindFlux( settings ); }
