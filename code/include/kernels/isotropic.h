#ifndef ISOTROPIC_H
#define ISOTROPIC_H

#include "scatteringkernelbase.h"

class Isotropic : public ScatteringKernel
{
  private:
    Isotropic() = delete;

  public:
    Isotropic( QuadratureBase* q );
    ~Isotropic();

    virtual Matrix GetScatteringKernel();
};

#endif
