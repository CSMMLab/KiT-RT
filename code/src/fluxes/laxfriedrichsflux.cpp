#include "fluxes/laxfriedrichsflux.h"

LaxFriedrichsFlux::LaxFriedrichsFlux( Config* settings ) : NumericalFlux( settings ) {}

double LaxFriedrichsFlux::Flux( const Vector& Omega, double psiL, double psiR, const Vector& n ) const {
    // double normN = norm( n );
    // return 0.5 * inner( n, Omega ) * ( psiL + psiR ) - 0.5 * ( psiR - psiL ) * norm( n ) / _dt;
    return -1.0;
}
