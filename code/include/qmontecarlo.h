#ifndef QMONTECARLO_H
#define QMONTECARLO_H

#include "quadrature.h"

class QMonteCarlo:public Quadrature
{
  public:
    QMonteCarlo( int order );
    virtual ~QMonteCarlo() {}


    std::string ComputeName();
    int ComputeNq();
    blaze::DynamicVector<blaze::DynamicVector<double>> ComputePoints();
    blaze::DynamicVector<double> ComputeWeights();
    blaze::DynamicVector<blaze::DynamicVector<int>> ComputeConnectivity();
};

#endif    // QMONTECARLO_H
