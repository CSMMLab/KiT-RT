#include "kernels/isotropic.h"

Isotropic::Isotropic( QuadratureBase* q ) : ScatteringKernel( q ) {}

Isotropic::~Isotropic() {}

SparseMatrix Isotropic::GetScatteringKernel() {
    unsigned nq = _q->GetNq();
    auto w      = _q->GetWeights();
    SparseMatrix kernel( nq, nq );
    for( unsigned i = 0; i < nq; ++i )
        for( unsigned j = 0; j < nq; ++j ) kernel( i, j ) = w[j] / ( 4 * M_PI );
    return kernel;
}
