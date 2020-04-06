#ifndef QLDFESA_H
#define QLDFESA_H

#include "qlookupquadrature.h"

class QLDFESA : public QLookupQuadrature
{
  public:
    QLDFESA( unsigned order );
    virtual ~QLDFESA() {}

    std::string ComputeName() override;
    VectorVectorU ComputeConnectivity() override;
};

#endif // QLDFESA_H
