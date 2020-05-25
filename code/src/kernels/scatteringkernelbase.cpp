#include "kernels/scatteringkernelbase.h"
#include "kernels/isotropic.h"

ScatteringKernel::ScatteringKernel( unsigned nq ) : _nq( nq ) {}

ScatteringKernel::~ScatteringKernel() {}

ScatteringKernel* ScatteringKernel::CreateScatteringKernel( KERNEL_NAME name, unsigned nq ) {
    switch( name ) {
        case KERNEL_Isotropic: return new Isotropic( nq );
        default: return new Isotropic( nq );
    }
}
