#ifndef SCATTERINGKERNELBASE_CPP
#define SCATTERINGKERNELBASE_CPP

#include "settings/globalconstants.h"
#include "typedef.h"

class ScatteringKernel
{
  private:
    ScatteringKernel() = delete;

  protected:
    unsigned _nq;

  public:
    ScatteringKernel( unsigned nq );
    ~ScatteringKernel();

    virtual SparseMatrix GetScatteringKernel() = 0;

    static ScatteringKernel* CreateScatteringKernel( KERNEL_NAME name, unsigned nq );
};

#endif    // SCATTERINGKERNELBASE_CPP
