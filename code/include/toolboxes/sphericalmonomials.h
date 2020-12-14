/*!
 * @file sphericalmonomials.h
 * @brief Class for efficient computation of a spherical monomial basis
 * @author S. Schotthöfer
 */

#ifndef SPHERICALMONOMIALS_H
#define SPHERICALMONOMIALS_H

#include "common/typedef.h"
#include "toolboxes/sphericalbase.h"
#include <vector>

class SphericalMonomials : public SphericalBase
{
  public:
    /*! @brief : Sets up class for monomial basis on sphere up to degree L.
     *           The basis then consists of N = L.
     *  @param : L_degree - maximum degree of spherical harmonics basis, 0 <= L <= 1000 (upper bound
     *                      due to numerical stability)
     * */
    SphericalMonomials( unsigned L_degree );

    /*! @brief  : Computes all N = L² +2L basis functions at point (my, phi)
     *  @param  : my = cos(theta) - spherical coordinate, -1 <= x <= 1
     *  @param  : phi - spherical coordinate, 0 <= phi <= 2*pi
     *  @return : vector of basis functions at point (my, phi) with size N = L² +2L
     */
    Vector ComputeSphericalBasis( double my, double phi ) override;

    /*! @brief  : Computes all N = L² +2L basis functions at point (x, y, z) on the unit sphere
     *  @param  : x,y,z = coordinates on unit sphere
     *  @return : vector of basis functions at point (x,y,z) with size N = L² +2L
     */
    Vector ComputeSphericalBasis( double x, double y, double z ) override;

    /*! @brief: Computes the length of the basis vector for a given max degree and
     *          spatial dimension dim. len of a single oder: (degree + _spatialDim -1) over (degree)
     *  @return: lenght of whole basis */
    unsigned GetBasisSize() override;

    /*! @brief: Computes the amount of lin. independent monomials of degree currDegreeL and
     *          spatial dimension dim. len of a single oder: (currDegreeL + _spatialDim -1) over (currDegreeL)
     *  @param: currDegreeL = degree of polynomials that are counted
     *  @return: lenght of a single dimension */
    unsigned GetCurrDegreeSize( unsigned currDegreeL ) override;

    /*! @brief: Computes global index of basis vector depending on order k and degree l
     *  @param: l = degree of polynomials l = 0,1,2,3,...
     *  @param: k = order of element of degree l. 0 <=k <=GetCurrDegreeSize(l) */
    unsigned GetGlobalIndexBasis( int l_degree, int k_order ) override;

  private:
    /*! @brief: Spatial dimension of the unit sphere (1,2,3) */
    unsigned _spatialDim;

    /*! @brief: spherical monomial basis function vector of
     *         degree  0 <= l <= L
     *         length : COmputed with ComputeBasisSize
     */
    Vector _YBasis;

    /*! @brief: Function to compute factorial of n (n!) in recursive manner */
    unsigned Factorial( unsigned n );

    /*! @brief: Function to compute first component of spherical unit vector
     *          Omega_x = sqrt(1-my*my)*cos(phi)
     *  @return: first component of spherical unit vector */
    double Omega_x( double my, double phi );

    /*! @brief: Function to compute first component of spherical unit vector
     *          Omega_x = sqrt(1-my*my)*sin(phi)
     *  @return: first component of spherical unit vector */
    double Omega_y( double my, double phi );

    /*! @brief: Function to compute first component of spherical unit vector
     *          Omega_z = my
     *  @return: first component of spherical unit vector */
    double Omega_z( double my );

    /*! @brief: Helper Function to compute basis^exponent. */
    double Power( double basis, unsigned exponent );
};
#endif    // SPHERICALMONOMIALS_H
