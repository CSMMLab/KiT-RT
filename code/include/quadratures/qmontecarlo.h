#ifndef QMONTECARLO_H
#define QMONTECARLO_H

#include "common/pch.h"

class Config;

class QMonteCarlo : public QuadratureBase
{
  public:
    QMonteCarlo( Config* settings );
    QMonteCarlo( unsigned quadOrder );
    inline ~QMonteCarlo() {}

    inline void SetName() override { _name = "Monte Carlo Quadrature."; }
    inline void SetNq() override { _nq = GetOrder(); }
    void SetPointsAndWeights() override;
    void SetConnectivity() override;
};

#endif    // QMONTECARLO_H
