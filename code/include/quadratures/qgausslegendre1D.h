#ifndef QGAUSSLEGENDRE1D_H
#define QGAUSSLEGENDRE1D_H

#include "quadraturebase.h"

class QGaussLegendre1D : public QuadratureBase
{
  private:
    double Pythag( const double a, const double b );
    std::pair<Vector, Matrix> ComputeEigenValTriDiagMatrix( const Matrix& mat );
    bool CheckOrder();

  public:
    QGaussLegendre1D( unsigned order );
    virtual ~QGaussLegendre1D() {}

    inline void SetName() override { _name = "Tensorized Gauss-Legendre quadrature."; }
    inline void SetNq() override { _nq = pow( GetOrder(), 2 ); }
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
};

#endif    // QGAUSSLEGENDRE1D_H
