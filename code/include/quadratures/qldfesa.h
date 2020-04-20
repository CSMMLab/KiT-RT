#ifndef QLDFESA_H
#define QLDFESA_H

#include "qlookupquadrature.h"

class QLDFESA : public QLookupQuadrature
{
  public:
    QLDFESA( unsigned order );
    virtual ~QLDFESA() {}

    inline void SetName() override { _name =  "LDFESA quadrature"; }
    void SetAvailOrders() override;
    void SetConnectivity() override;

    std::string GetLookupTable() override;
};

#endif // QLDFESA_H
