#include "snsolver.h"
#include <mpi.h>

SNSolver::SNSolver( Settings* settings ) : Solver( settings ) {}

void SNSolver::Solve() {

    // TODO: set up initial condition

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _psi;

    // loop over energies (pseudo-time)
    for( unsigned n = 0; n < _nTimeSteps; ++n ) {
        std::cout << "time " << n * _dt << std::endl;
        // loop over all spatial cells
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] ) continue;
            // loop over all ordinates
            for( unsigned k = 0; k < _nq; ++k ) {
                psiNew[j][k] = 0.0;
                // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                for( unsigned l = 0; l < _neighbors[j].size(); ++l ) {
                    // store flux contribution on psiNew_sigmaSH20 to save memory
                    psiNew[j][k] -= _g->Flux( _quadPoints[k], _psi[j][k], _psi[_neighbors[j][l]][k], _normals[j][l] );
                    std::cout << _g->Flux( _quadPoints[k], _psi[j][k], _psi[_neighbors[j][l]][k], _normals[j][l] ) << " "
                              << _g->Flux( _quadPoints[k], _psi[_neighbors[j][l]][k], _psi[j][k], -_normals[j][l] ) << std::endl;
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[j][k] = _psi[j][k] + ( _dt / _areas[j] ) * psiNew[j][k] + _dt * _sigmaTH20[n] * _psi[j][k];
            }
            // compute scattering effects
            psiNew[j] += _sigmaSH20[n] * _psi[j] * _weights;    // multiply scattering matrix with psi
        }
        _psi = psiNew;
    }
}

void SNSolver::SolveMPI() {
    /*
        // setup MPI variables: mype is PE index, npes is number of PEs, ierr is MPI message flag
        int mype, npes, ierr;

        // determine mype and npes
        ierr = MPI_Comm_rank( MPI_COMM_WORLD, &mype );
        ierr = MPI_Comm_size( MPI_COMM_WORLD, &npes );

        // determine size of quadrature array
        int nqPE = int( ( _nq - 1 ) / npes ) + 1;
        // nqPEMax needed for allocation
        int nqPEMax = nqPE;
        if( mype == npes - 1 ) {
            nqPE = _nq - mype * nqPE;
            if( nqPE < 0 ) {
                nqPE = 0;
            }
        }

        int master = 0;    // define PE 0 to be master
        int tag    = 0;
        int kStart = mype * ( ( _nq - 1 ) / npes + 1.0 );
        int kEnd   = kStart + nqPE - 1;

    // TODO: store local psi: dim(psi) = (NCells,nqPE)
    // TODO: store local scattering matrix: dim(_sigmaSH20) = (nTimeSteps,_nq,nqPE)
    // TODO: store local weights: dim(_weights) = (nqPE)

    // TODO: set up initial condition

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _psi;

    // loop over energies (pseudo-time)
    for( unsigned n = 0; n < _nTimeSteps; ++n ) {
        // loop over all spatial cells
        for( unsigned j = 0; j < _NCells; ++j ) {

            // loop over all ordinates
            for( unsigned k = 0; k < unsigned( nqPE ); ++k ) {
                psiNew[j][k] = 0.0;
                // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                for( unsigned l = 0; l < _neighbors[j].size(); ++l ) {
                    // store flux contribution on psiNew to save memory
                    psiNew[j][k] -= ( _dt / _areas[j] ) * _g->Flux( _quadPoints[k], _psi[j][k], _psi[_neighbors[j][l]][k], _normals[j][l] );
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[j][k] = _psi[j][k] + psiNew[j][k] + _dt * _sigmaTH20[n] * _psi[j][k];
            }
            // compute scattering effects
            Vector scattering = _sigmaSH20[n] * _psi[j] * _weights;    // multiply scattering matrix with dim(scattering) = _nq -> improve later
            psiNew[j] += scattering;
        }
        _psi = psiNew;
        // psiNew.reset();
    }*/
}

void SNSolver::Save() const {
    std::vector<std::string> fieldNames{"flux"};
    std::vector<double> flux( _nCells, 0.0 );
    for( unsigned i = 0; i < _nCells; ++i ) {
        flux[i] = dot( _psi[i], _weights );
    }
    std::vector<std::vector<double>> scalarField( 1, flux );
    std::vector<std::vector<std::vector<double>>> results{scalarField};
    ExportVTK( _settings->GetOutputFile(), results, fieldNames, _settings, _mesh );
}
