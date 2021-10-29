/*!
 * @file isotropic1D.h
 * @brief Class for computing an isotropic scattering kernel of kinetic equations
 * @author ?
 */

#ifndef ISOTROPIC_H
#define ISOTROPIC_H

#include "scatteringkernelbase.hpp"

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
