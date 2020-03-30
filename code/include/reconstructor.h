#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include "settings.h"
#include "typedef.h"

class Reconstructor
{
  public:
    /**
     * @brief Reconstruction
     * @param settings
     */
    Reconstructor( Settings* settings );

    /**
     * @brief Flux computes flux on edge for fixed ordinate at a given edge
     * @param Omega fixed ordinate for flux computation
     * @param psiL left solution state
     * @param psiR right solution state
     * @param n scaled normal vector of given edge
     * @return reconstructed slope
     */
    virtual double Slope( const Vector& Omega, double psiL, double psiR, const Vector& n ) const;
};

#endif    // RECONSTRUCTOR_H
