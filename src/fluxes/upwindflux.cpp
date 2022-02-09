#include "fluxes/upwindflux.hpp"

UpwindFlux::UpwindFlux() : NumericalFlux() {}

double UpwindFlux::Flux( const Vector& Omega, double psiL, double psiR, const Vector& n ) const {
    double inner = Omega[0] * n[0] + Omega[1] * n[1];    // Only use x and y axis in 2d case
    if( inner > 0 ) {
        return inner * psiL;
    }
    else {
        return inner * psiR;
    }
}

Vector UpwindFlux::Flux( const Matrix AxPlus,
                         const Matrix AxMinus,
                         const Matrix AyPlus,
                         const Matrix AyMinus,
                         const Matrix /*AzPlus*/,
                         const Matrix /*AzMinus*/,
                         const Vector psiL,
                         const Vector psiR,
                         const Vector n ) const {
    // 2d only atm!!!

    Vector resultFlux( psiR.size(), 0 );
    // x dir
    if( n[0] > 0 ) {
        resultFlux += n[0] * AxPlus * psiL + n[0] * AxMinus * psiR;
    }
    else {
        resultFlux += n[0] * AxPlus * psiR + n[0] * AxMinus * psiL;
    }
    // y dir
    if( n[1] > 0 ) {
        resultFlux += n[1] * AyPlus * psiL + n[1] * AyMinus * psiR;
    }
    else {
        resultFlux += n[1] * AyPlus * psiR + n[1] * AyMinus * psiL;
    }

    return resultFlux;
}

Vector UpwindFlux::FluxXZ( const Matrix AxPlus,
                           const Matrix AxMinus,
                           const Matrix /*AyPlus*/,
                           const Matrix /*AyMinus*/,
                           const Matrix AzPlus,
                           const Matrix AzMinus,
                           const Vector psiL,
                           const Vector psiR,
                           const Vector n ) const {
    // 2d only atm!!!

    Vector resultFlux( psiR.size(), 0 );
    // x dir
    if( n[0] > 0 ) {
        resultFlux += n[0] * AxPlus * psiL + n[0] * AxMinus * psiR;
    }
    else {
        resultFlux += n[0] * AxPlus * psiR + n[0] * AxMinus * psiL;
    }
    // y dir
    if( n[1] > 0 ) {
        resultFlux += n[1] * AzPlus * psiL + n[1] * AzMinus * psiR;
    }
    else {
        resultFlux += n[1] * AzPlus * psiR + n[1] * AzMinus * psiL;
    }

    return resultFlux;
}

void UpwindFlux::FluxVanLeer( const Matrix& Ax,
                              const Matrix& AxAbs,
                              const Matrix& /*Ay*/,
                              const Matrix& /*AyAbs*/,
                              const Matrix& Az,
                              const Matrix& AzAbs,
                              const Vector& psiL,
                              const Vector& psiR,
                              const Vector& n,
                              Vector& resultFlux ) const {

    // resultFlux += n[0] * ( AxPlus * psiL + AxMinus * psiR ) + n[1] * ( AyPlus * psiL + AyMinus * psiR );
    resultFlux += 0.5 * ( n[0] * ( Ax * ( psiL + psiR ) - AxAbs * ( psiR - psiL ) ) + n[1] * ( Az * ( psiL + psiR ) - AzAbs * ( psiR - psiL ) ) );
}
