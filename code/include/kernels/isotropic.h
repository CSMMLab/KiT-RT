#ifndef ISOTROPIC_H
#define ISOTROPIC_H

#include "scatteringkernelbase.h"

class Isotropic : public ScatteringKernel
{
  private:
    Isotropic() = delete;

  public:
    Isotropic( unsigned nq );
    ~Isotropic();

    virtual SparseMatrix GetScatteringKernel();
};

#endif
