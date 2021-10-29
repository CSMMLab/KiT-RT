#include "kernels/scatteringkernelbase.hpp"
#include "kernels/isotropic.hpp"
#include "kernels/isotropic1D.hpp"

ScatteringKernel::ScatteringKernel( QuadratureBase* quad ) : _quad( quad ) {}

ScatteringKernel::~ScatteringKernel() {}

ScatteringKernel* ScatteringKernel::CreateScatteringKernel( KERNEL_NAME name, QuadratureBase* quad ) {
    switch( name ) {
        case KERNEL_Isotropic: return new Isotropic( quad );
        case KERNEL_Isotropic1D: return new Isotropic( quad );
        default: return new Isotropic( quad );
    }
}
