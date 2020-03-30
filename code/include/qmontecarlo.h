#ifndef QMONTECARLO_H
#define QMONTECARLO_H

#include "quadrature.h"

class QMonteCarlo : public Quadrature
{
  public:
    QMonteCarlo( unsigned order );
    virtual ~QMonteCarlo() {}

    std::string ComputeName();
    unsigned ComputeNq();
    VectorVector ComputePoints();
    Vector ComputeWeights();
    VectorVectorU ComputeConnectivity();
};

#endif    // QMONTECARLO_H
