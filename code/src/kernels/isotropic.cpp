#include "kernels/isotropic.h"

Isotropic::Isotropic( unsigned nq ) : ScatteringKernel( nq ) {}

Isotropic::~Isotropic() {}

SparseMatrix Isotropic::GetScatteringKernel() {
    SparseMatrix kernel( _nq, _nq );
    for( unsigned i = 0; i < _nq; ++i )
        for( unsigned j = 0; j < _nq; ++j ) kernel( i, j ) = 1 / ( 4 * M_PI );
    return kernel;
}
