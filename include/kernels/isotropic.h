/*!
 * @file isotropic1D.h
 * @brief Class for computing an isotropic scattering kernel of kinetic equations
 * @author ?
 */

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
    virtual Matrix GetScatteringKernelFirstCollision( unsigned nqF, Vector weightsF );
};

#endif
