/*!
 * @file isotropic1D.h
 * @brief Class for computing an isotropic scattering kernel of kinetic equations
 * @author ?
 */

#ifndef ISOTROPIC_H
#define ISOTROPIC_H

#include "common/typedef.h"    // for Matrix
#include "scatteringkernelbase.h"

class QuadratureBase;

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
