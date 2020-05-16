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
    virtual double Flux( const Vector& Omega, double psiL, double psiR, const Vector& n ) const;
};

#endif    // UPWINDFLUX_H
