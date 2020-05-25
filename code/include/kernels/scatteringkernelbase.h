#ifndef SCATTERINGKERNELBASE_CPP
#define SCATTERINGKERNELBASE_CPP

#include "quadratures/quadraturebase.h"
#include "settings/globalconstants.h"
#include "typedef.h"

class ScatteringKernel
{
  private:
    ScatteringKernel() = delete;

  protected:
    QuadratureBase* _q;

  public:
    ScatteringKernel( QuadratureBase* q );
    ~ScatteringKernel();

    virtual SparseMatrix GetScatteringKernel() = 0;

    static ScatteringKernel* CreateScatteringKernel( KERNEL_NAME name, QuadratureBase* q );
};

#endif    // SCATTERINGKERNELBASE_CPP
