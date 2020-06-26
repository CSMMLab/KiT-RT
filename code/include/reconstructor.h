#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include "settings/config.h"
#include "settings/typedef.h"

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
    
};

#endif    // RECONSTRUCTOR_H


double FortSign( double a, double b ) ;
double LMinMod( double sL, double sR ) ;
double LVanLeer( double sL, double sR ) ;
double LSuperBee( double sL, double sR ) ;
double LVanAlbaba( double sL, double sR ) ;
double LWENOJS( double x ) ;