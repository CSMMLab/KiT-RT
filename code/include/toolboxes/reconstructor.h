/*!
 * @file: reconstructor.h
 * @brief Class to create second order (in space) schemes for the advection solver.
 *         This class is currently unused. But in the future, the second order capabilities of the code
 *         may be stored here.
 * @author: T. Xiao
 */

#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <math.h>    // for fabs, fmin, fmax
#include <string>

class Config;

class Reconstructor
{
  protected:
    unsigned _reconsOrder;

  public:
    /**
     * @brief Reconstruction
     * @param settings
     */
    Reconstructor( Config* settings );
    static Reconstructor* Create( Config* settings );

    unsigned inline GetReconsOrder() { return _reconsOrder; }

    /*! Method 1: structured developing
     * @brief Slope of angular flux psi inside a given cell
     * @param uL left solution state
     * @param uC center solution state
     * @param uR right solution state
     * @param dxL left slope
     * @param dxR right slope
     * @param limiter name of the applied limiter (invalid argument results in disabled limiter)
     * @return reconstructed slope
     */
    virtual double ReconstructSlopeStruct( double uL, double uC, double uR, double dxL, double dxR, std::string limiter ) const;
};

#endif    // RECONSTRUCTOR_H

double FortSign( double a, double b );
double LMinMod( double sL, double sR );
double LVanLeer( double sL, double sR );
double LSuperBee( double sL, double sR );
double LVanAlbaba( double sL, double sR );
double LWENOJS( double x );
