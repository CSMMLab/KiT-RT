#ifndef QRECTANGULAR
#define QRECTANGULAR

#include "quadraturebase.hpp"

class QRectangular : public QuadratureBase
{
    // Implementation is done accordingly to Kendall Atkinson 1981, Australian Matematical Society.
  public:
    QRectangular( Config* settings );
    virtual ~QRectangular() {}

protected:
    inline void SetName() override { _name = "Rectangular Midpoint quadrature"; }
    void SetNq() override;
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
    void ScalePointsAndWeights(double velocityScaling) override;
};

class QRectangular2D : public QuadratureBase
{
  public:
    QRectangular2D( Config* settings );
    virtual ~QRectangular2D() {}

protected:
    inline void SetName() override { _name = "Rectangular Midpoint quadrature 2D projection"; }
    void SetNq() override;
    void CheckOrder();
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
    void ScalePointsAndWeights(double velocityScaling) override;
};

class QRectangular1D : public QuadratureBase
{
  public:
    QRectangular1D( Config* settings );
    virtual ~QRectangular1D() {}
  protected:
    inline void SetName() override { _name = "Rectangular quadrature 1D projection"; }
    void SetNq() override;
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
    void ScalePointsAndWeights(double velocityScaling) override;
};
#endif    // QRECTANGULAR
