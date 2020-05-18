#ifndef LAXFRIEDRICHSFLUX_H
#define LAXFRIEDRICHSFLUX_H

#include "numericalflux.h"
#include "typedef.h"

class LaxFriedrichsFlux : public NumericalFlux
{
    double _dt;

  public:
    /**
     * @brief LaxFriedrichsFlux
     * @param settings
     */
    LaxFriedrichsFlux( Config* settings );

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

#endif    // LAXFRIEDRICHSFLUX_H
