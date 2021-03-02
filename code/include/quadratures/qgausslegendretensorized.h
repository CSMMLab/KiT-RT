#ifndef QGAUSSLEGENDRETENSORIZED_H
#define QGAUSSLEGENDRETENSORIZED_H

#include "common/pch.h"

class Config;

class QGaussLegendreTensorized : public QuadratureBase
{
    // Implementation is done accordingly to Kendall Atkinson 1981, Australian Matematical Society.
  private:
    double Pythag( const double a, const double b );
    std::pair<Vector, Matrix> ComputeEigenValTriDiagMatrix( const Matrix& mat );
    bool CheckOrder();

  public:
    QGaussLegendreTensorized( Config* settings );

    virtual ~QGaussLegendreTensorized() {}

    inline void SetName() override { _name = "Tensorized Gauss-Legendre quadrature"; }
    void SetNq() override;
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
};

#endif    // QGAUSSLEGENDRETENSORIZED_H
