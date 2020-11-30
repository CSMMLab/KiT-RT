#ifndef PRODUCTQUADRATURE_H
#define PRODUCTQUADRATURE_H

#include "quadraturebase.h"

class ProductQuadrature : public QuadratureBase
{
    // Implementation is done accordingly to Kendall Atkinson 1981, Australian Matematical Society.
  private:
    double Pythag( const double a, const double b );
    std::pair<Vector, Matrix> ComputeEigenValTriDiagMatrix( const Matrix& mat );
    bool CheckOrder();

  public:
    ProductQuadrature( Config* settings );
    ProductQuadrature( unsigned order );

    virtual ~ProductQuadrature() {}

    inline void SetName() override { _name = "Product quadrature"; }
    inline void SetNq() override { _nq = 4 * pow( GetOrder(), 2 ); }
    void SetPointsAndWeights() override;
    void SetConnectivity() override;

    inline VectorVector GetPointsSphere() const override { return _pointsSphere; }
};

#endif    // PRODUCTQUADRATURE_H
