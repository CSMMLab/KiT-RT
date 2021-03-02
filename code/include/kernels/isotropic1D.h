/*!
 * @file isotropic1D.h
 * @brief Class for computing an isotropic scattering kernel of 1D kinetic equations
 * @author ?
 */

#ifndef ISOTROPIC1D_H
#define ISOTROPIC1D_H

#include "common/pch.h"

class QuadratureBase;

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
