#include "numericalflux.h"
#include "upwindflux.h"

NumericalFlux::NumericalFlux( CConfig* settings ) {}

NumericalFlux* NumericalFlux::Create( CConfig* settings ) { return new UpwindFlux( settings ); }
