#include "fluxes/laxfriedrichsflux.h"
#include "fluxes/numericalflux.h"    // for NumericalFlux

LaxFriedrichsFlux::LaxFriedrichsFlux() : NumericalFlux() {}

double LaxFriedrichsFlux::Flux( const Vector& /*Omega*/, double /*psiL*/, double /*psiR*/, const Vector& /*n*/ ) const {
    // double normN = norm( n );
    // return 0.5 * inner( n, Omega ) * ( psiL + psiR ) - 0.5 * ( psiR - psiL ) * norm( n ) / _dt;
    return -1.0;
}
