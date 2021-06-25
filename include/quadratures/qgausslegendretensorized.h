#ifndef QGAUSSLEGENDRETENSORIZED_H
#define QGAUSSLEGENDRETENSORIZED_H

#include "quadraturebase.h"

class QGaussLegendreTensorized : public QuadratureBase
{
    // Implementation is done accordingly to Kendall Atkinson 1981, Australian Matematical Society.
  protected:
    double Pythag( const double a, const double b );
    std::pair<Vector, Matrix> ComputeEigenValTriDiagMatrix( const Matrix& mat );
    virtual bool CheckOrder();

  public:
    QGaussLegendreTensorized( Config* settings );

    virtual ~QGaussLegendreTensorized() {}

    virtual inline void SetName() override { _name = "Tensorized Gauss-Legendre quadrature"; }
    void SetNq() override;
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
};

class QGaussLegendreTensorized2D : public QGaussLegendreTensorized
{
    // Implementation is done accordingly to Kendall Atkinson 1981, Australian Matematical Society.
  private:
    bool CheckOrder() override;

  public:
    QGaussLegendreTensorized2D( Config* settings );

    virtual ~QGaussLegendreTensorized2D() {}

    inline void SetName() override { _name = "Tensorized Gauss-Legendre quadrature 2D projection"; }
    void SetNq() override;
    void SetPointsAndWeights() override;
};
#endif    // QGAUSSLEGENDRETENSORIZED_H
