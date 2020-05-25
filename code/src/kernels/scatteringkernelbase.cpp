#include "kernels/scatteringkernelbase.h"
#include "kernels/isotropic.h"

ScatteringKernel::ScatteringKernel( QuadratureBase* q ) : _q( q ) {}

ScatteringKernel::~ScatteringKernel() {}

ScatteringKernel* ScatteringKernel::CreateScatteringKernel( KERNEL_NAME name, QuadratureBase* q ) {
    switch( name ) {
        case KERNEL_Isotropic: return new Isotropic( q );
        default: return new Isotropic( q );
    }
}
