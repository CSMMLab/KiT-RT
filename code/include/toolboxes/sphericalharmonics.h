/*!
 * @file sphericalharmonics.h
 * @brief Class for efficient computation of a spherical harmonics basis
 *        according to the algorithm described in
 *        "Associated Legendre Polynomials and Spherical Harmonics Computation  for Chemistry Applications"
 *        by Limpanupab, T. and Milthorpe, J. (https://arxiv.org/pdf/1410.1748.pdf)
 * @author S. Schotthöfer
 *
 */

#ifndef SPHERICALHARMONICS_H
#define SPHERICALHARMONICS_H

#include "common/typedef.h"
#include "toolboxes/sphericalbase.h"
#include <vector>

class SphericalHarmonics : public SphericalBase
{
  public:
    /*! @brief Sets up class for spherical harmonics basis based on legendre
     *         polynoms and associated legendre polynoms up to degree L.
     *         The basis then consists of N = L² +2L basis functions.
     *  @param L_degree - maximum degree of spherical harmonics basis, 0 <= L <= 1000 (upper bound
     *                      due to numerical stability)
     * */
    SphericalHarmonics( unsigned L_degree );

    /*! @brief  Computes all N = L² +2L basis functions at point (my, phi)
     *  @param  my = cos(theta) - spherical coordinate, -1 <= x <= 1
     *  @param  phi - spherical coordinate, 0 <= phi <= 2*pi
     *  @return vector of basis functions at point (my, phi) with size N = L² +2L
     */
    Vector ComputeSphericalBasis( double my, double phi ) override;

    /*! @brief   Computes all N = L² +2L basis functions at point (x, y, z) on the unit sphere
     *  @param   x coordinates on unit sphere
     *  @param   y coordinates on unit sphere
     *  @param   z coordinates on unit sphere
     *  @return  vector of basis functions at point (x,y,z) with size N = L² +2L
     */
    Vector ComputeSphericalBasis( double x, double y, double z ) override;

    /*! @brief  Computes an entire set of (komplex congjugate) P_l^k and stores
     *          it in the vector _assLegendreP
     *  @param  my = cos(theta)  - spherical coordinate, -1 <=  my <= 1
     *  @return  Associated Legendre Polynomial at my for all l and k
     */
    std::vector<double> GetAssLegendrePoly( const double my );

    /*! @brief Returns length of the basis, i.e. number of elements of the basis */
    unsigned GetBasisSize() override;

    /*! @brief Return number of basis functions with degree equals to currDegree
     *  @param currDegree degree of polynomials that are counted   */
    unsigned GetCurrDegreeSize( unsigned currDegree ) override;

    /*! @brief  helper function to get the global index for given k and l of
     *          the basis function Y_k^l.
     *  @param  l_degree - current degree of basis function, 0 <= l <= L
     *  @param  k_order  - current order of basis function,  -l <= k <= l      */
    unsigned GetGlobalIndexBasis( int l_degree, int k_order ) override;

  private:
    /*! @brief maximal degree of the spherical harmonics basis (this is "L" in the comments)*/
    unsigned _LMaxDegree;

    /*! @brief coefficients for the computations of the basis
     *         length of _aParam, _bParam : L + (L*(L+1))/2
     */
    std::vector<double> _aParam;
    std::vector<double> _bParam;

    /*! @brief komplex conjugate of the associated legendre polynomials of
     *         degree  0 <= l <= L and order 0 <= k <= l
     *        length of _assLegendreP : L + (L*(L+1))/2
     */
    std::vector<double> _assLegendreP;

    /*! @brief spherical harmonic basis functions of
     *         degree  0 <= l <= L and order -l <= k <= l
     *         length : N = L² +2L
     */
    Vector _YBasis;

    /*! @brief  helper function to get the global index for given k and l of
     *          the associated legendre polynomial P_k^l.
     *  @param  l_degree - current degree of basis function, 0 <= l <= L
     *  @param  k_order  - current order of basis function,  0 <= k <= l
     */
    unsigned inline GlobalIdxAssLegendreP( unsigned l_degree, unsigned k_order ) { return k_order + ( l_degree * ( l_degree + 1 ) ) / 2; }

    /*! @brief  computes values of a_param and b_param
     */
    void ComputeCoefficients();

    /*! @brief  Computes an entire set of (komplex congjugate) P_l^k and stores
     *          it in the vector _assLegendreP
     *  @param  my  equals cos(theta)  - spherical coordinate, -1 <=  my <= 1
     */
    void ComputeAssLegendrePoly( const double my );

    /*! @brief Computes the spherical harmonics basis function up to degree _LmaxDegree at
     *          polar coordinates (theta, psi) and stores the result in _YBasis;
     *  @param phi spherical coordinate
     */
    void ComputeYBasis( const double phi );
};
#endif    // SPHERICALHARMONICS_H
