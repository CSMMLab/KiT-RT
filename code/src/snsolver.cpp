#include "snsolver.h"

SNSolver::SNSolver( Settings* settings ) : Solver( settings ) {}

void SNSolver::Solve() {
    // TODO main SN time loop

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    Matrix psiNew = _psi;

    // loop over energies (pseudo-time)
    for( unsigned n = 0; n < _nTimeSteps; ++n ) {
        // loop over all spatial cells
        for( unsigned j = 0; j < _NCells; ++j ) {
            // loop over all ordinates
            for( unsigned k = 0; k < _nq; ++k ) {
                // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                for( unsigned l = 0; l < _neighbors[j].size(); ++l ) {
                    psiNew( j, k ) -= ( _dt / _areas[j] ) * _g->Flux( _quadPoints[k], _psi( j, k ), _psi( _neighbors[j][l], k ), _normals[j][l] );
                }
                psiNew( j, k ) = _psi( j, k ) + psiNew( j, k );
            }
        }
        _psi = psiNew;
        psiNew.reset();
    }
}
