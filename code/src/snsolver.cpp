#include "snsolver.h"

SNSolver::SNSolver( Config* settings ) : Solver( settings ) {}

void SNSolver::Solve() {
    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _psi;
    double dFlux        = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "dFlux" );

    // loop over energies (pseudo-time)
    // for( unsigned n = 0; n < _nEnergies; ++n ) {
    unsigned ctr = 0;
    while( dFlux > 1e-6 ) {
        // loop over all spatial cells
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            // loop over all ordinates
            for( unsigned k = 0; k < _nq; ++k ) {
                psiNew[j][k] = 0.0;
                // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                for( unsigned l = 0; l < _neighbors[j].size(); ++l ) {
                    // store flux contribution on psiNew_sigmaS to save memory
                    psiNew[j][k] -= _g->Flux( _quadPoints[k], _psi[j][k], _psi[_neighbors[j][l]][k], _normals[j][l] );
                }
                // time update angular flux with numerical flux and total scattering cross section
                // psiNew[j][k] = _psi[j][k] + ( _dE / _areas[j] ) * psiNew[j][k] - _dE * _sigmaT[n][j] * _psi[j][k];
                psiNew[j][k] = _psi[j][k] + ( _dE / _areas[j] ) * psiNew[j][k] - _dE * _sigmaT[0][j] * _psi[j][k];
            }
            // compute scattering effects
            // psiNew[j] += _dE * _scatteringKernel * _sigmaS[n][j] * _psi[j] * _weights;    // multiply scattering matrix with psi
            psiNew[j] += _dE * _scatteringKernel * _sigmaS[0][j] * _psi[j] * _weights;    // multiply scattering matrix with psi
            // add external source contribution
            // psiNew[j] += _Q[n][j] * _areas[j];
            psiNew[j] += _Q[0][j] * _areas[j];
        }
        _psi = psiNew;
        for( unsigned i = 0; i < _nCells; ++i ) {
            fluxNew[i] = dot( _psi[i], _weights );
        }
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        // if( rank == 0 ) log->info( "{:03.8f}   {:01.5e}", _energies[n], dFlux );
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e}", ++ctr * _dE, dFlux );
    }
}

void SNSolver::Save() const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<double> flux( _nCells, 0.0 );
    for( unsigned i = 0; i < _nCells; ++i ) {
        flux[i] = dot( _psi[i], _weights );
    }
    std::vector<std::vector<double>> scalarField( 1, flux );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile(), results, fieldNames, _mesh );
    auto log = spdlog::get( "event" );
    log->info( "Result successfully exported to '{0}'!", _settings->GetOutputFile() );
}
