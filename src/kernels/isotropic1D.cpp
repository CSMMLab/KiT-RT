#include "kernels/isotropic1D.h"
#include "quadratures/quadraturebase.h"
#include "toolboxes/errormessages.h"

Isotropic1D::Isotropic1D( QuadratureBase* quad ) : ScatteringKernel( quad ) {}

Isotropic1D::~Isotropic1D() {}

Matrix Isotropic1D::GetScatteringKernel() {
    unsigned nq = _quad->GetNq();
    auto w      = _quad->GetWeights();

    Matrix kernel( nq, nq );
    for( unsigned i = 0; i < nq; ++i )
        for( unsigned j = 0; j < nq; ++j ) kernel( i, j ) = w[j] / 2;

    // scale kernel to ensure mass conservation
    double tmp;
    for( unsigned i = 0; i < nq; ++i ) {
        tmp = 0.0;
        for( unsigned j = 0; j < nq; ++j ) {
            tmp += kernel( i, j );
        }
        for( unsigned j = 0; j < nq; ++j ) {
            // kernel( i, j ) /= tmp;
        }
    }
    return kernel;
}

Matrix Isotropic1D::GetScatteringKernelFirstCollision( unsigned , Vector ){
    ErrorMessages::Error("First Collision ScatteringKernel not defined for Isotropic1D case. Please choose Isotropic Kernel!", CURRENT_FUNCTION );
} // not used
