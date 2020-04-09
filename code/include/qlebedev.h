#ifndef QLEBEDEV_H
#define QLEBEDEV_H

#include "qlookupquadrature.h"

class QLebedev : public QLookupQuadrature
{
  public:
    QLebedev( unsigned order );
    virtual ~QLebedev() {}

    inline void SetName() override { _name =  "Lebedev quadrature"; }
    void SetAvailOrders() override;
    void SetDataInfo() override;
    void SetConnectivity() override;

};

#endif // QLEBEDEV_H
