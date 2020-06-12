#include "fluxes/upwindflux.h"

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
    // std::cout << "AxPlus" << AxPlus << std::endl;
    // std::cout << "psiL" << psiL << std::endl;
    // Vector temp = AxPlus * psiL;
    // std::cout << "AxPlus * psiL" << temp << std::endl;
    // std::cout << "n_x *AxPlus * psiL" << n[0] * temp << std::endl;

    resultFlux += ( n[0] * AxPlus + n[1] * AyPlus ) * psiL + ( n[0] * AxMinus + n[1] * AyMinus ) * psiR;
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
                              const Matrix& Ay,
                              const Matrix& AyAbs,
                              const Matrix& Az,
                              const Matrix& AzAbs,
                              const Vector& psiL,
                              const Vector& psiR,
                              const Vector& n,
                              Vector& resultFlux ) const {

    // resultFlux += n[0] * ( AxPlus * psiL + AxMinus * psiR ) + n[1] * ( AyPlus * psiL + AyMinus * psiR );
    resultFlux += 0.5 * ( n[0] * ( Ax * ( psiL + psiR ) - AxAbs * ( psiR - psiL ) ) + n[1] * ( Az * ( psiL + psiR ) - AzAbs * ( psiR - psiL ) ) );
}
