#ifndef QLEVELSYMMETRIC_H
#define QLEVELSYMMETRIC_H

#include "qlookupquadrature.h"

class QLevelSymmetric : public QLookupQuadrature
{
  public:
    QLevelSymmetric( unsigned order );
    virtual ~QLevelSymmetric() {}

    std::string ComputeName() override;
    VectorVectorU ComputeConnectivity() override;
};

#endif    // QLEVELSYMMETRIC_H
