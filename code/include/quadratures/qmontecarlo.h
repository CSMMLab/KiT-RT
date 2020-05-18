#ifndef QMONTECARLO_H
#define QMONTECARLO_H

#include "quadraturebase.h"

class QMonteCarlo : public QuadratureBase
{
  public:
    QMonteCarlo( unsigned order );
    virtual ~QMonteCarlo() {}

    inline void SetName() override { _name = "Monte Carlo Quadrature."; }
    inline void SetNq() override { _nq = GetOrder(); }
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
};

#endif    // QMONTECARLO_H
