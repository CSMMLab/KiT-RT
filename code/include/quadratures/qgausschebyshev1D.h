/*!
 * @file config.h
 * @brief Class to compute 1D Gauss Legendre Quadrature
 * @author J. Kusch
 *
 * Disclaimer: This class structure was copied and modifed with open source permission from SU2 v7.0.3 https://su2code.github.io/
 */

#ifndef QGAUSSCHEBYSHEV_H
#define QGAUSSCHEBYSHEV_H

#include "quadraturebase.h"

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
     *  @param double(f)( double x0, double x1, double x2 ) : density function that depends on a three spatial dimensions.
     *  @returns double result: result of the quadrature rule */
    double Integrate( double( f )( double x0, double x1, double x2 ) ) override;

    /*! @brief Integrates f(x,y,z) with the quadrature.
     *  @param double(f)( double my, double phi ) : density function that depends on a spherical coordinates.
     *  @returns double result: result of the quadrature rule */
    double IntegrateSpherical( double( f )( double my, double phi ) ) override;

    /*! @brief Integrates vector valued f(x,y,z) with the quadrature. Each dimension is integrated by itself.
     *  @param : double(f)( double x0, double x1, double x2 ) : density function that depends on a three spatial dimensions.
     *  @param :  len : lenght of vector
     *  @returns double result: result of the quadrature rule (vector valued) */
    std::vector<double> Integrate( std::vector<double>( f )( double x0, double x1, double x2 ), unsigned /* len */ ) override;
};

#endif    // QGAUSSCHEBYSHEV_H
