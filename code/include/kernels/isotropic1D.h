#ifndef ISOTROPIC1D_H
#define ISOTROPIC1D_H

#include "scatteringkernelbase.h"

class Isotropic1D : public ScatteringKernel
{
  private:
    Isotropic1D() = delete;

  public:
    Isotropic1D( QuadratureBase* q );
    ~Isotropic1D();

    virtual Matrix GetScatteringKernel();
};

#endif
