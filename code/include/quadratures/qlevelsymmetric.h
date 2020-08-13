#ifndef QLEVELSYMMETRIC_H
#define QLEVELSYMMETRIC_H

#include "qlookupquadrature.h"

class QLevelSymmetric : public QLookupQuadrature
{
  public:
    QLevelSymmetric( Config* settings );
    virtual ~QLevelSymmetric() {}

    inline void SetName() override { _name = "Level Symmetric quadrature"; }
    void SetAvailOrders() override;
    void SetConnectivity() override;

    std::string GetLookupTable() override;
};

#endif    // QLEVELSYMMETRIC_H
