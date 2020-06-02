#ifndef UPWINDFLUX_H
#define UPWINDFLUX_H

#include "numericalflux.h"
#include "typedef.h"

class UpwindFlux : public NumericalFlux
{
  public:
    /**
     * @brief UpwindFlux
     * @param settings
     */
    UpwindFlux( Config* settings );

    /**
     * @brief Flux computes flux on edge for fixed ordinate at a given edge
     * @param Omega fixed ordinate for flux computation
     * @param psiL left solution state
     * @param psiR right solution state
     * @param n scaled normal vector of given edge
     * @return numerical flux value
     */
    double Flux( const Vector& Omega, double psiL, double psiR, const Vector& n ) const override;

    /**
     * @brief Flux      : Computes <Steger Warming> upwinding scheme for given flux jacobians of the PN Solver at a given edge and stores it in
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
    void Flux( const Matrix& AxPlus,
               const Matrix& AxMinus,
               const Matrix& AyPlus,
               const Matrix& AyMinus,
               const Matrix& AzPlus,
               const Matrix& AzMinus,
               const Vector& psiL,
               const Vector& psiR,
               const Vector& n,
               Vector& resultFlux ) const override;
};

#endif    // UPWINDFLUX_H
