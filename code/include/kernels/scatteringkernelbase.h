/*!
 * @file scatteringkernelbase.h
 * @brief Base class for computing the scattering kernel of kinetic equations
 * @author ?
 */

#ifndef SCATTERINGKERNELBASE_CPP
#define SCATTERINGKERNELBASE_CPP

#include "common/globalconstants.h"
#include "common/typedef.h"

// Forward Declaration
class QuadratureBase;

class ScatteringKernel
{
  private:
    ScatteringKernel() = delete;

  protected:
    QuadratureBase* _quad; /*! @brief: Pointer to the quadrature used to compute the scattering integral */

  public:
    /*! @brief: Copies the pointer of quad. Does not create an own quad */
    ScatteringKernel( QuadratureBase* quad );

    virtual ~ScatteringKernel();

    /*! @brief: Computes the scattering kernel and for the whole SN system and stores it in a Matrix
        @return: Matrix with discretized scattering kernel */
    virtual Matrix GetScatteringKernel() = 0;

    /*! @brief: Creates an object of the child class of ScatteringKernelBase corresponding to the enum KERNEL_NAME */
    static ScatteringKernel* CreateScatteringKernel( KERNEL_NAME name, QuadratureBase* quad );
};

#endif    // SCATTERINGKERNELBASE_CPP
