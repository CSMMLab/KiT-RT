/*!
 * @file config.h
 * @brief Class to compute 1D Gauss Legendre Quadrature
 * @author J. Kusch
 *
 * Disclaimer: This class structure was copied and modifed with open source permission from SU2 v7.0.3 https://su2code.github.io/
 */

#ifndef QGAUSSCHEBYSHEV_H
#define QGAUSSCHEBYSHEV_H

#include "quadraturebase.hpp"
#include <iostream>

class QGaussChebyshev1D : public QuadratureBase
{
  private:
    bool CheckOrder();

  public:
    QGaussChebyshev1D( Config* settings );
    QGaussChebyshev1D( unsigned quadOrder );
    virtual ~QGaussChebyshev1D() {}

    inline void SetName() override { _name = "Gauss-Chebychev quadrature 1D."; }
    inline void SetNq() override { _nq = _order; }
    void SetPointsAndWeights() override;
    void SetConnectivity() override;

    /*! @brief Integrates f(x,y,z) with the quadrature.
     *  @param f density function that depends on a three spatial dimensions.
     *  @returns result of the quadrature rule */
    double Integrate( double ( *f )( double, double, double ) ) override;

    /*! @brief Integrates f(x,y,z) with the quadrature.
     *  @param f density function that depends on a spherical coordinates.
     *  @returns result of the quadrature rule */
    double IntegrateSpherical( double ( *f )( double, double ) ) override;

    /*! @brief Integrates vector valued f(x,y,z) with the quadrature. Each dimension is integrated by itself.
     *  @param f density function that depends on a three spatial dimensions.
     *  @returns result of the quadrature rule (vector valued) */
    std::vector<double> Integrate( std::vector<double> ( *f )( double, double, double ), unsigned /* len */ ) override;
};

#endif    // QGAUSSCHEBYSHEV_H
