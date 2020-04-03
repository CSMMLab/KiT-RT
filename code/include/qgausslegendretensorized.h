#ifndef QGAUSSLEGENDRETENSORIZED_H
#define QGAUSSLEGENDRETENSORIZED_H

#include "quadrature.h"

class QGaussLegendreTensorized : public Quadrature
{
  public:
    QGaussLegendreTensorized( unsigned order );
    virtual ~QGaussLegendreTensorized() {}

    std::string ComputeName() override;
    unsigned ComputeNq() override;
    VectorVector ComputePoints() override;
    Vector ComputeWeights() override;
    VectorVectorU ComputeConnectivity() override;
};

#endif    // QGAUSSLEGENDRETENSORIZED_H
