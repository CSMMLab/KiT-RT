#ifndef QLEBEDEV_H
#define QLEBEDEV_H

#include "qlookupquadrature.h"

class QLebedev : public QLookupQuadrature
{
  public:
    QLebedev( Config* settings );
    virtual ~QLebedev() {}

    inline void SetName() override { _name = "Lebedev quadrature"; }
    void SetAvailOrders() override;
    void SetConnectivity() override;

    std::string GetLookupTable() override;
};

#endif    // QLEBEDEV_H
