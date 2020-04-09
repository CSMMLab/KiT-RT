#ifndef QMONTECARLO_H
#define QMONTECARLO_H

#include "quadrature.h"

class QMonteCarlo : public Quadrature
{
  public:
    QMonteCarlo( unsigned order );
    virtual ~QMonteCarlo() {}

    inline void SetName() override { _name =  "Monte Carlo Quadrature."; }
    void SetNq() override;
    void SetPointsAndWeights() override;
    void SetConnectivity()override;
};

#endif    // QMONTECARLO_H
