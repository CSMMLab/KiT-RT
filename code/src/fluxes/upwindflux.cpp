#include "fluxes/upwindflux.h"

UpwindFlux::UpwindFlux() : NumericalFlux() {}

double UpwindFlux::Flux( const Vector& Omega, double psiL, double psiR, const Vector& n ) const {
    double inner = Omega[0] * n[0] + Omega[1] * n[1];
    if( inner > 0 ) {
        return inner * psiL;
    }
    else {
        return inner * psiR;
    }
}

/**
 * @brief Flux      : Computes <Linear> upwinding scheme for given flux jacobians of the PN Solver at a given edge and stores it in
 *                    resultFlux
 * @param AxPlus    : Positive part of the flux jacobian in x direction
 * @param AxMinus   : Negative part of the flux jacobian in x direction
 * @param AyPlus    : Positive part of the flux jacobian in y direction
 * @param AyMinus   : Negative part of the flux jacobian in y direction
 * @param AzPlus    : Positive part of the flux jacobian in z direction
 * @param AzMinus   : Negative part of the flux jacobian in z direction
 * @param psiL      : Solution state of left hand side control volume
 * @param psiR      : Solution state of right hand side control volume
 * @param n         : Normal vector at the edge between left and right control volume
 * @return resultFlux: Vector with resulting flux.
 */

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

/**
 * @brief Flux      : Computes <VanLeer> upwinding scheme for given flux jacobians of the PN Solver at a given edge and stores it in
 *                    resultFlux
 * @param AxPlus    : Positive part of the flux jacobian in x direction
 * @param AxMinus   : Negative part of the flux jacobian in x direction
 * @param AyPlus    : Positive part of the flux jacobian in y direction
 * @param AyMinus   : Negative part of the flux jacobian in y direction
 * @param AzPlus    : Positive part of the flux jacobian in z direction
 * @param AzMinus   : Negative part of the flux jacobian in z direction
 * @param psiL      : Solution state of left hand side control volume
 * @param psiR      : Solution state of right hand side control volume
 * @param n         : Normal vector at the edge between left and right control volume
 * @param resultFlux: Vector with resulting flux.
 * @return          : void
 */
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
