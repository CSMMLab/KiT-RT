#ifndef QMONTECARLO_H
#define QMONTECARLO_H

#include "quadraturebase.h"

class QMonteCarlo : public QuadratureBase
{
  public:
    QMonteCarlo( unsigned order );
    inline ~QMonteCarlo() {}

    inline void SetName() override { _name = "Monte Carlo Quadrature."; }
    inline void SetNq() override { _nq = GetOrder(); }
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
    VectorVector GetPointsSphere() const override;
};

#endif    // QMONTECARLO_H
