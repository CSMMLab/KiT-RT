/*!
 * @file sphericalmonomials.h
 * @brief Class for efficient computation of a spherical monomial basis
 * @author S. Schotthöfer
 */

#ifndef SPHERICALMONOMIALS_H
#define SPHERICALMONOMIALS_H

#include "common/typedef.hpp"
#include "velocitybasis/sphericalbase.hpp"
#include <vector>

class SphericalMonomials : public SphericalBase
{
  public:
    /*! @brief Sets up class for monomial basis on sphere up to degree L.
     *           The basis then consists of N = L.
     *  @param L_degree maximum degree of spherical harmonics basis, 0 <= L <= 1000 (upper bound
     *                      due to numerical stability)
     * */
    SphericalMonomials( unsigned L_degree );

    /*! @brief Sets up class for monomial basis on sphere up to degree L.
     *         The basis then consists of N = L.
     *  @param L_degree maximum degree of spherical harmonics basis, 0 <= L <= 1000 (upper bound
     *                    due to numerical stability)
     *  @param spatialDim spatial dimensioniality of the simulation
     * */
    SphericalMonomials( unsigned L_degree, unsigned short spatialDim );

    /*! @brief  Computes all N = L² +2L basis functions at point (my, phi)
     *  @param  my cos(theta) - spherical coordinate, -1 <= x <= 1
     *  @param  phi spherical coordinate, 0 <= phi <= 2*pi
     *  @param  r radius of sphere, default is 1, i.e. unit sphere
     *  @return vector of basis functions at point (my, phi) with size N = L² +2L
     */
    Vector ComputeSphericalBasis( double my, double phi, double r = 1.0 ) override;

    /*! @brief Computes all N = L² +2L basis functions at point (x, y, z) on the unit sphere
     *  @param x = coordinates on sphere
     *  @param y = coordinates on sphere
     *  @param z = coordinates on sphere
     *  @return vector of basis functions at point (x,y,z) with size N = L² +2L
     */
    Vector ComputeSphericalBasisKarthesian( double x, double y, double z ) override;

    /*! @brief Computes the length of the basis vector for a given max degree and
     *          spatial dimension dim. len of a single oder: (degree + _spatialDim -1) over (degree)
     *  @return lenght of whole basis */
    unsigned GetBasisSize() override;

    /*! @brief Computes the amount of lin. independent monomials of degree currDegreeL and
     *          spatial dimension dim. len of a single oder: (currDegreeL + _spatialDim -1) over (currDegreeL)
     *  @param currDegreeL = degree of polynomials that are counted
     *  @return lenght of a single dimension */
    unsigned GetCurrDegreeSize( unsigned currDegreeL ) override;

    /*! @brief Computes global index of basis vector depending on order k and degree l
     *  @param l_degree degree of polynomials l = 0,1,2,3,...
     *  @param k_order order of element of degree l. 0 <=k <=GetCurrDegreeSize(l) */
    unsigned GetGlobalIndexBasis( int l_degree, int k_order ) override;

  private:
    /*! @brief spherical monomial basis function vector of
     *         degree  0 <= l <= L
     *         length : COmputed with ComputeBasisSize
     */
    Vector _YBasis;

    /*! @brief Function to compute factorial of n (n!) in recursive manner */
    unsigned Factorial( unsigned n );

    /*! @brief Function to compute first component of spherical unit vector
     *          Omega_x = r* sqrt(1-my*my)*cos(phi)
     *  @return first component of spherical unit vector */
    double Omega_x( double my, double phi, double r = 1.0 );

    /*! @brief Function to compute first component of spherical unit vector
     *          Omega_x = r* sqrt(1-my*my)*sin(phi)
     *  @return first component of spherical unit vector */
    double Omega_y( double my, double phi, double r = 1.0 );

    /*! @brief Function to compute first component of spherical unit vector
     *          Omega_z = r* my
     *  @return first component of spherical unit vector */
    double Omega_z( double my, double r = 1.0 );

    /*! @brief Helper Function to compute basis^exponent. */
    double Power( double basis, unsigned exponent );

    /*! @brief Helper to compute the basis in 1D, 2D or 3D, depending on choice of _spatialDim */
    Vector ComputeSphericalBasis1D( double my, double r = 1.0 );             /*!< @brief Only Z achsis of the 3D case */
    Vector ComputeSphericalBasis2D( double my, double phi, double r = 1.0 ); /*!< @brief Only X and Y achsis of the 3D case */
    Vector ComputeSphericalBasis3D( double my, double phi, double r = 1.0 );
};
#endif    // SPHERICALMONOMIALS_H
