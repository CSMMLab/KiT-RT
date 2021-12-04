#ifndef QMIDPOINTTENSORIZED
#define QMIDPOINTTENSORIZED

#include "quadraturebase.hpp"

class QMidpointTensorized : public QuadratureBase
{
    // Implementation is done accordingly to Kendall Atkinson 1981, Australian Matematical Society.
  public:
    QMidpointTensorized( Config* settings );
    virtual ~QMidpointTensorized() {}

protected:
    inline void SetName() override { _name = "Tensorized Midpoint quadrature"; }
    void SetNq() override;
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
};

class QMidpointTensorized2D : public QuadratureBase
{
  public:
    QMidpointTensorized2D( Config* settings );
    virtual ~QMidpointTensorized2D() {}

protected:
    inline void SetName() override { _name = "Tensorized Midpoint quadrature 2D projection"; }
    void SetNq() override;
    void CheckOrder();
    void SetPointsAndWeights() override;
    void SetConnectivity() override;

};

class QMidpoint1D : public QuadratureBase
{
  public:
    QMidpoint1D( Config* settings );
    virtual ~QMidpoint1D() {}
  protected:
    inline void SetName() override { _name = "Midpoint quadrature 1D projection"; }
    void SetNq() override;
    void SetPointsAndWeights() override;
    void SetConnectivity() override;

};
#endif    // QMIDPOINTTENSORIZED
