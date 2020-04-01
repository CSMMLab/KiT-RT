#include "upwindflux.h"

UpwindFlux::UpwindFlux( Settings* settings ) : NumericalFlux( settings ) {}

double UpwindFlux::Flux( const Vector& Omega, double psiL, double psiR, const Vector& n ) const {
    if( inner( Omega, n ) > 0 ) {
        return inner( Omega, n ) * psiL;
    }
    else {
        return inner( Omega, n ) * psiR;
    }
}
