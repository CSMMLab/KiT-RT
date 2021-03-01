#include "kernels/isotropic.h"

#include "blaze/math/dense/DynamicVector.h"    // for DynamicVector
#include "kernels/scatteringkernelbase.h"      // for ScatteringKernel
#include "quadratures/quadraturebase.h"
#include <math.h>    // for M_PI

Isotropic::Isotropic( QuadratureBase* quad ) : ScatteringKernel( quad ) {}

Isotropic::~Isotropic() {}

Matrix Isotropic::GetScatteringKernel() {
    unsigned nq = _quad->GetNq();
    auto w      = _quad->GetWeights();
    Matrix kernel( nq, nq );
    for( unsigned i = 0; i < nq; ++i )
        for( unsigned j = 0; j < nq; ++j ) kernel( i, j ) = w[j] / ( 4 * M_PI );

    // scale kernel to ensure mass conservation
    double tmp;
    for( unsigned i = 0; i < nq; ++i ) {
        tmp = 0.0;
        for( unsigned j = 0; j < nq; ++j ) {
            tmp += kernel( i, j );
        }
        for( unsigned j = 0; j < nq; ++j ) {
            kernel( i, j ) /= tmp;
        }
    }
    return kernel;
}
