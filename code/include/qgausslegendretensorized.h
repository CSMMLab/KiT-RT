#ifndef QGAUSSLEGENDRETENSORIZED_H
#define QGAUSSLEGENDRETENSORIZED_H

#include "quadrature.h"

class QGaussLegendreTensorized : public Quadrature
{
  public:
    QGaussLegendreTensorized( unsigned order );
    virtual ~QGaussLegendreTensorized() {}

    std::string ComputeName();
    unsigned ComputeNq();
    VectorVector ComputePoints();
    Vector ComputeWeights();
    VectorVectorU ComputeConnectivity();
};

#endif    // QGAUSSLEGENDRETENSORIZED_H
