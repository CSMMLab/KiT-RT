
#include "pnsolver.h"
#include <mpi.h>

PNSolver::PNSolver( Config* settings ) : Solver( settings ) {}

void PNSolver::Solve() {
    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _psi;

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "#placeholder" );

    // loop over energies (pseudo-time)
    for( unsigned n = 0; n < _nTimeSteps; ++n ) {
        // loop over all spatial cells
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] ) continue;
            // loop over all moments
            for( unsigned k = 0; k < _nq; ++k ) {
                psiNew[j][k] = 0.0;
                // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                for( unsigned l = 0; l < _neighbors[j].size(); ++l ) {
                    // store flux contribution on psiNew_sigmaSH20 to save memory
                    // TODO A matrix!
                    psiNew[j][k] -= _g->Flux( _quadPoints[k], _psi[j][k], _psi[_neighbors[j][l]][k], _normals[j][l] );
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[j][k] = _psi[j][k] + ( _dt / _areas[j] ) * psiNew[j][k] + _dt * _sigmaTH20[n] * _psi[j][k];
            }
            // compute scattering effects
            psiNew[j] += _sigmaSH20[n] * _psi[j] * _weights;    // multiply scattering matrix with psi
        }
        _psi = psiNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e}", n * _dt, 0.0 );
    }
}
