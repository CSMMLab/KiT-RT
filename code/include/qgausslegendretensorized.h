#ifndef QGAUSSLEGENDRETENSORIZED_H
#define QGAUSSLEGENDRETENSORIZED_H

#include "quadrature.h"

class QGaussLegendreTensorized:public Quadrature
{
  public:
    QGaussLegendreTensorized( int order );
    virtual ~QGaussLegendreTensorized() {}


    std::string ComputeName();
    int ComputeNq();
    blaze::DynamicVector<blaze::DynamicVector<double>> ComputePoints();
    blaze::DynamicVector<double> ComputeWeights();
    blaze::DynamicVector<blaze::DynamicVector<int>> ComputeConnectivity();
};

#endif // QGAUSSLEGENDRETENSORIZED_H

