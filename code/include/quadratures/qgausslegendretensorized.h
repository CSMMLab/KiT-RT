#ifndef QGAUSSLEGENDRETENSORIZED_H
#define QGAUSSLEGENDRETENSORIZED_H

#include "quadraturebase.h"

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

    inline VectorVector GetPointsSphere() const override { return _pointsSphere; }
};

#endif    // QGAUSSLEGENDRETENSORIZED_H
