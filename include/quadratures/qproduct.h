/*! @file: qproduct.h
 *  @brief Product quadrature implementation. Implementation is done accordingly to Kendall Atkinson 1981, Australian Matematical Society.
 *  @author: J. Kusch
 */

#ifndef PRODUCTQUADRATURE_H
#define PRODUCTQUADRATURE_H

#include "quadraturebase.h"

class QProduct : public QuadratureBase
{
  private:
    double Pythag( const double a, const double b );
    std::pair<Vector, Matrix> ComputeEigenValTriDiagMatrix( const Matrix& mat );
    bool CheckOrder();

  public:
    QProduct( Config* settings );
    QProduct( unsigned order );

    virtual ~QProduct() {}

    inline void SetName() override { _name = "Product quadrature"; }
    inline void SetNq() override { _nq = 4 * pow( GetOrder(), 2 ); }
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
};

#endif    // PRODUCTQUADRATURE_H
