#ifndef NUMERICALFLUX_H
#define NUMERICALFLUX_H

#include "typedef.h"

class NumericalFlux
{
public:
    NumericalFlux(Settings* settings);
    static NumericalFlux* Create( Settings* settings );
    /**
     * @brief Flux computes flux on edge for fixed ordinate at a given edge
     * @param Omega fixed ordinate for flux computation
     * @param psiL left solution state
     * @param psiR right solution state
     * @param n scaled normal vector of given edge
     * @return numerical flux value
     */
    virtual double Flux(const Vector& Omega, double psiL, double psiR, const Vector& n)const = 0;
};

#endif // NUMERICALFLUX_H
