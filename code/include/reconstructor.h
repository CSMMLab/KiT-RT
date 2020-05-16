#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include "settings/config.h"
#include "typedef.h"

class Reconstructor
{
  public:
    /**
     * @brief Reconstruction
     * @param settings
     */
    Reconstructor( Config* settings );

    /** Method 1: structured developing
     * @brief Slope of angular flux psi inside a given cell
     * @param Omega fixed ordinate for flux computation
     * @param psiL left solution state
     * @param psiR right solution state
     * @param n scaled normal vector of given edge
     * @return reconstructed slope
     */

    /** Method 2: unstructured developing
     * @brief Slope of angular flux psi inside a given cell
     * @param Omega fixed ordinate for flux computation
     * @param psiL left solution state
     * @param psiR right solution state
     * @param n scaled normal vector of given edge
     * @return reconstructed slope
     */

    virtual double ReconstructSlopeStruct( double uL, double uC, double uR, double dxL, double dxR, std::string limiter ) const;
    virtual void ReconstructSlopeUnstruct( unsigned nq, unsigned ncell, VectorVector& psi, VectorVector& psiDerX, VectorVector& psiDerY, Vector& area, VectorVector& neighbor, VectorVector& nx, VectorVector& ny ) const;

};

#endif // RECONSTRUCTOR_H
