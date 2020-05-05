#include "upwindflux.h"

UpwindFlux::UpwindFlux( CConfig* settings ) : NumericalFlux( settings ) {}

double UpwindFlux::Flux( const Vector& Omega, double psiL, double psiR, const Vector& n ) const {
    // if( inner( Omega, n ) > 0 ) {
    double inner = Omega[0] * n[0] + Omega[1] * n[1];
    if( inner > 0 ) {
        return inner * psiL;
    }
    else {
        return inner * psiR;
    }
}
