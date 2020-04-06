#ifndef QLEBEDEV_H
#define QLEBEDEV_H

#include "qlookupquadrature.h"

class QLebedev : public QLookupQuadrature
{
  public:
    QLebedev( unsigned order );
    virtual ~QLebedev() {}

    std::string ComputeName() override;
    VectorVectorU ComputeConnectivity() override;
};

#endif // QLEBEDEV_H
