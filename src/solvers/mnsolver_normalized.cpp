#include "solvers/mnsolver_normalized.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "entropies/entropybase.hpp"
#include "fluxes/numericalflux.hpp"
#include "optimizers/optimizerbase.hpp"
#include "problems/problembase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "solvers/mnsolver.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/sphericalbase.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

// externals
#include "spdlog/spdlog.h"
#include <iostream>
#include <mpi.h>

MNSolverNormalized::MNSolverNormalized( Config* settings ) : MNSolver( settings ) { _u0 = Vector( _nCells, 0.0 ); }

MNSolverNormalized::~MNSolverNormalized() {}

void MNSolverNormalized::IterPreprocessing( unsigned /*idx_pseudotime*/ ) {

    // ------- Entropy closure Step ----------------
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        _u0[idx_cell] = _sol[idx_cell][0];
        _sol[idx_cell] /= _u0[idx_cell];    // assume _u0 > 0 always!!
    }
    TextProcessingToolbox::PrintVectorVectorToFile( _sol, "solution", _nCells, _nSystem );

    // TextProcessingToolbox::PrintVectorVector( _sol );
    _optimizer->SolveMultiCell( _alpha, _sol, _momentBasis );

    // ------- Solution reconstruction step ----
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            // compute the kinetic density at all grid cells
            _kineticDensity[idx_cell][idx_quad] =
                _u0[idx_cell] * _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], _momentBasis[idx_quad] ) );
        }
        _sol[idx_cell] *= _u0[idx_cell];
        if( _settings->GetRealizabilityReconstruction() ) ComputeRealizableSolution( idx_cell );
    }

    // ------ Compute slope limiters and cell gradients ---
    if( _reconsOrder > 1 ) {
        _mesh->ComputeSlopes( _nq, _solDx, _solDy, _kineticDensity );               // parallel
        _mesh->ComputeLimiter( _nq, _solDx, _solDy, _kineticDensity, _limiter );    // parallel
    }
}
