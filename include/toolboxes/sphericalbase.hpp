/*!
 * @file sphericalbase.h
 * @brief Base Class to handle basis classes on the unit sphere
 * @author S. Schotthöfer
 *
 */

#ifndef SPHERICALBASE_H
#define SPHERICALBASE_H

#include "common/typedef.hpp"
class Config;

class SphericalBase
{
  public:
    SphericalBase() {}
    virtual ~SphericalBase() {}

    /*! @brief Create a set of basis functions on the unit sphere defined in settings
     *  @param settings pointer to the config object
     *  @returns Pointer to the createt basis class
     */
    static SphericalBase* Create( Config* settings );

    /*! @brief  Computes all N basis functions at point (my, phi)
     *  @param  my cos(theta) - spherical coordinate, -1 <= x <= 1
     *  @param  phi spherical coordinate, 0 <= phi <= 2*pi
     *  @param  r radius of sphere. Hardcoded to r=1 for spherical harmonics.
     *  @return vector of basis functions at point (my, phi) with size N
     */
    virtual Vector ComputeSphericalBasis( double my, double phi, double r = 1.0 ) = 0;

    /*! @brief  Computes all basis functions at point (x, y, z) on the unit sphere
     *  @param  x,y,z coordinates on unit sphere
     *  @return vector of basis functions at point (x,y,z) with size N
     */
    virtual Vector ComputeSphericalBasisKarthesian( double x, double y, double z ) = 0;

    /*! @brief Return size of complete Basisvector */
    virtual unsigned GetBasisSize() = 0;

    /*! @brief Return number of basis functions with degree equals to currDegree
     *  @param currDegree must be smaller equals _LMaxDegree     */
    virtual unsigned GetCurrDegreeSize( unsigned currDegree ) = 0;

    /*! @brief Computes global index of basis vector depending on order k and degree l
     *  @param l_degree = degree of polynomials l = 0,1,2,3,...
     *  @param k_order = order of element of degree l. !ATTENTION. Requirements are different for monomials and harmonics! */
    virtual unsigned GetGlobalIndexBasis( int l_degree, int k_order ) = 0;

  protected:
    /*! @brief maximal (polynomial) degree of the spherical basis (this is "L" in the comments)*/
    unsigned _LMaxDegree;
    /*! @brief Spatial dimension of the unit sphere (projection) (1,2,3) */
    unsigned _spatialDim;
};

#endif    // SPHERICALBASE_H
