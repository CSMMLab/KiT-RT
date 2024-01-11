#include "solvers/mnsolver_normalized.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "entropies/entropybase.hpp"
#include "fluxes/numericalflux.hpp"
#include "optimizers/newtonoptimizer.hpp"
#include "optimizers/optimizerbase.hpp"
#include "optimizers/partregularizednewtonoptimizer.hpp"
#include "problems/problembase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "solvers/mnsolver.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "velocitybasis/sphericalbase.hpp"

// externals
#include "spdlog/spdlog.h"
#include <iostream>

MNSolverNormalized::MNSolverNormalized( Config* settings ) : MNSolver( settings ) {
    _u0 = Vector( _nCells, 0.0 );
    if( _settings->GetRegularizerGamma() > 0 ) {
        _optimizer2 = new PartRegularizedNewtonOptimizer( settings );
    }
    else {
        _optimizer2 = new NewtonOptimizer( settings );
    }
    _sol2 = VectorVector( _nCells, Vector( _nSystem, 0.0 ) );
}

MNSolverNormalized::~MNSolverNormalized() {}

void MNSolverNormalized::IterPreprocessing( unsigned idx_pseudotime ) {
    Vector alpha_norm_per_cell( _nCells, 0 );    // ONLY FOR DEBUGGING! THIS SLOWS DOWN THE CODE

    // if (idx_pseudotime < 1600){
    //     // preprocessor
    //     _optimizer3->SolveMultiCell( _alpha, _sol, _momentBasis, alpha_norm_per_cell );    // Newton for the first few iterations
    //
    //    #pragma omp parallel for
    //    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
    //        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
    //            // compute the kinetic density at all grid cells
    //            _kineticDensity[idx_cell][idx_quad] = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], _momentBasis[idx_quad] ) );
    //        }
    //    }
    //}
    // else
    {
        // std::cout << "Start Network based closure";
        // if (idx_pseudotime==1600){
        //     #pragma omp parallel for
        //     for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        //          _sol[idx_cell][0] = max(1e-10, _sol[idx_cell][0]);
        //     }
        // }
        //  ------- Entropy closure Step ----------------
#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
            _u0[idx_cell] = _sol[idx_cell][0];
            _sol[idx_cell] /= _u0[idx_cell];    // assume _u0 > 0 always!!
        }

        _optimizer->SolveMultiCell( _alpha, _sol, _momentBasis, alpha_norm_per_cell );    // Newton for the first few iterations

        int norm_err = 0.0;
        // Check if solution is close to reconstruction, if not, apply newton optimizer
#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
            _optimizer2->ReconstructMoments( _sol2[idx_cell], _alpha[idx_cell], _momentBasis );
            //_optimizer2->ReconstructMoments( _sol2[idx_cell], _alpha[idx_cell], _momentBasis );
            if( blaze::norm( _sol2[idx_cell] - _sol[idx_cell] ) / blaze::norm( _sol[idx_cell] ) > _settings->GetNewtonOptimizerEpsilon() ) {
                //_optimizer2->Solve( _alpha[idx_cell], _sol[idx_cell], _momentBasis, idx_cell );    // Newton postprocessor
                norm_err++;
            }
            // Dynamic closure
            alpha_norm_per_cell[idx_cell] = 0.0;
            for( unsigned idx_sys = 1; idx_sys < _nSystem; idx_sys++ ) {
                alpha_norm_per_cell[idx_cell] += _alpha[idx_cell][idx_sys] * _alpha[idx_cell][idx_sys];
            }
        }
        std::cout << "Number of inaccurate cells: " << norm_err << "\n";

        // ------- Solution reconstruction step ----
#pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
            alpha_norm_per_cell[idx_cell] *= _momentBasis[0][0] * 0.5 * _settings->GetRegularizerGamma();    // is constant
            if( _settings->GetEntropyDynamicAnsatz() ) {
                for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                    // compute the kinetic density at all grid cells
                    _kineticDensity[idx_cell][idx_quad] =
                        _u0[idx_cell] *
                        _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], _momentBasis[idx_quad] ) - alpha_norm_per_cell[idx_cell] );
                }
            }
            else {
                for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                    // compute the kinetic density at all grid cells
                    _kineticDensity[idx_cell][idx_quad] =
                        _u0[idx_cell] * _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], _momentBasis[idx_quad] ) );
                }
            }
            if( _settings->GetRealizabilityReconstruction() ) {
                ComputeRealizableSolution( idx_cell );
            }
            _sol[idx_cell] *= _u0[idx_cell];
        }
    }

    // ------ Compute slope limiters and cell gradients ---
    if( _reconsOrder > 1 ) {
        _mesh->ComputeSlopes( _nq, _solDx, _solDy, _kineticDensity );               // parallel
        _mesh->ComputeLimiter( _nq, _solDx, _solDy, _kineticDensity, _limiter );    // parallel
    }
}
