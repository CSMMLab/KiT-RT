/*!
 * @file sphericalbase.h
 * @brief Base Class to handle basis classes on the unit sphere
 * @author S. Schotthöfer
 *
 */

#ifndef SPHERICALBASE_H
#define SPHERICALBASE_H

#include "common/typedef.h"
class Config;

class SphericalBase
{
  public:
    SphericalBase() {}
    ~SphericalBase() {}

    /*! @brief: Create a set of basis functions on the unit sphere defined in settings
     *  @param: Pointer to the config file
     *  @returns: Pointer to the createt basis class */
    static SphericalBase* Create( Config* settings );

    /*! @brief  : Computes all N basis functions at point (my, phi)
     *  @param  : my = cos(theta) - spherical coordinate, -1 <= x <= 1
     *  @param  : phi - spherical coordinate, 0 <= phi <= 2*pi
     *  @return : vector of basis functions at point (my, phi) with size N
     */
    virtual Vector ComputeSphericalBasis( double my, double phi ) = 0;

    /*! @brief  : Computes all basis functions at point (x, y, z) on the unit sphere
     *  @param  : x,y,z = coordinates on unit sphere
     *  @return : vector of basis functions at point (x,y,z) with size N
     */
    virtual Vector ComputeSphericalBasis( double x, double y, double z ) = 0;

    virtual unsigned GetBasisSize() = 0;
};

#endif    // SPHERICALBASE_H
