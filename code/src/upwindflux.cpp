#include "upwindflux.h"

UpwindFlux::UpwindFlux( Config* settings ) : NumericalFlux( settings ) {}

double UpwindFlux::Flux( const Vector& Omega, double psiL, double psiR, const Vector& n ) const {
    double inner = Omega[0] * n[0] + Omega[1] * n[1];
    if( inner > 0 ) {
        return inner * psiL;
    }
    else {
        return inner * psiR;
    }
}

void UpwindFlux::Flux( const Matrix& AxPlus,
                       const Matrix& AxMinus,
                       const Matrix& AyPlus,
                       const Matrix& AyMinus,
                       const Matrix& AzPlus,
                       const Matrix& AzMinus,
                       const Vector& psiL,
                       const Vector& psiR,
                       const Vector& n,
                       Vector& resultFlux ) const {
    // 2d only atm!!!
    resultFlux = ( n[0] * AxPlus + n[1] * AyPlus ) * psiL + ( n[0] * AxMinus + n[1] * AyMinus ) * psiR;
}
