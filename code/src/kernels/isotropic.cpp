#include "kernels/isotropic.h"

Isotropic::Isotropic( QuadratureBase* quad ) : ScatteringKernel( quad ) {}

Isotropic::~Isotropic() {}

Matrix Isotropic::GetScatteringKernel() {
    unsigned nq = _quad->GetNq();
    auto w      = _quad->GetWeights();
    Matrix kernel( nq, nq );
    for( unsigned i = 0; i < nq; ++i )
        for( unsigned j = 0; j < nq; ++j ) kernel( i, j ) = w[j] / ( 4 * M_PI );
    return kernel;
}
