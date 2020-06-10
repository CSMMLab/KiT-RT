#ifndef SCATTERINGKERNELBASE_CPP
#define SCATTERINGKERNELBASE_CPP

#include "quadratures/quadraturebase.h"
#include "settings/globalconstants.h"
#include "settings/typedef.h"

class ScatteringKernel
{
  private:
    ScatteringKernel() = delete;

  protected:
    QuadratureBase* _quad;

  public:
    ScatteringKernel( QuadratureBase* quad );
    ~ScatteringKernel();

    virtual Matrix GetScatteringKernel() = 0;

    static ScatteringKernel* CreateScatteringKernel( KERNEL_NAME name, QuadratureBase* quad );
};

#endif    // SCATTERINGKERNELBASE_CPP
