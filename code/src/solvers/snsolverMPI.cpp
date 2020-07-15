#include "solvers/snsolverMPI.h"
#include "common/config.h"
#include "common/io.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

SNSolverMPI::SNSolverMPI( Config* settings ) : SNSolver( settings ) {}

void SNSolverMPI::Solve() {
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
    // TODO: store local scattering matrix: dim(_sigmaS) = (nTimeSteps,_nq,nqPE)
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
                psiNew[j][k] = _psi[j][k] + psiNew[j][k] + _dt * _sigmaT[n] * _psi[j][k];
            }
            // compute scattering effects
            Vector scattering = _sigmaS[n] * _psi[j] * _weights;    // multiply scattering matrix with dim(scattering) = _nq -> improve later
            psiNew[j] += scattering;
        }
        _psi = psiNew;
        // psiNew.reset();
    }*/
}

void SNSolverMPI::Save() const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<double> flux( _nCells, 0.0 );
    for( unsigned i = 0; i < _nCells; ++i ) {
        flux[i] = dot( _sol[i], _weights );
    }
    std::vector<std::vector<double>> scalarField( 1, flux );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile(), results, fieldNames, _mesh );
    auto log = spdlog::get( "event" );
    log->info( "Result successfully exported to '{0}'!", _settings->GetOutputFile() );
}
