#ifndef QMONTECARLO_H
#define QMONTECARLO_H

#include "quadrature.h"

class QMonteCarlo : public Quadrature
{
  public:
    QMonteCarlo( unsigned order );
    virtual ~QMonteCarlo() {}

    std::string ComputeName() override;
    unsigned ComputeNq() override;
    VectorVector ComputePoints() override;
    Vector ComputeWeights() override;
    VectorVectorU ComputeConnectivity() override;
};

#endif    // QMONTECARLO_H
