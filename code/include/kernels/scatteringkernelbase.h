#ifndef SCATTERINGKERNELBASE_CPP
#define SCATTERINGKERNELBASE_CPP

#include "settings/globalconstants.h"
#include "settings/typedef.h"

// Forward Declaration
class QuadratureBase;

class ScatteringKernel
{
  private:
    ScatteringKernel() = delete;

  protected:
    QuadratureBase* _quad;

  public:
    ScatteringKernel( QuadratureBase* quad );
    virtual ~ScatteringKernel();

    virtual Matrix GetScatteringKernel() = 0;

    static ScatteringKernel* CreateScatteringKernel( KERNEL_NAME name, QuadratureBase* quad );
};

#endif    // SCATTERINGKERNELBASE_CPP
