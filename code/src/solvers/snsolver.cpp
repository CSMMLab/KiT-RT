#include "solvers/snsolver.h"

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
    for( unsigned idx_energy = 0; idx_energy < _nEnergies; ++idx_energy ) {
        // loop over all spatial cells
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
            // loop over all ordinates
            for( unsigned idx_ord = 0; idx_ord < _nq; ++idx_ord ) {
                psiNew[idx_cell][idx_ord] = 0.0;
                // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); ++idx_neighbor ) {
                    // store flux contribution on psiNew_sigmaS to save memory
                    if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                        continue;    // adiabatic wall, add nothing
                    // psiNew[idx_cell][idx_ord] +=
                    //    _g->Flux( _quadPoints[idx_ord], _psi[idx_cell][idx_ord], _psi[idx_cell][idx_ord], _normals[idx_cell][idx_neighbor] );
                    else
                        psiNew[idx_cell][idx_ord] += _g->Flux( _quadPoints[idx_ord],
                                                               _psi[idx_cell][idx_ord],
                                                               _psi[_neighbors[idx_cell][idx_neighbor]][idx_ord],
                                                               _normals[idx_cell][idx_neighbor] );
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[idx_cell][idx_ord] = _psi[idx_cell][idx_ord] - ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_ord] -
                                            _dE * _sigmaT[idx_energy][idx_cell] * _psi[idx_cell][idx_ord];
            }
            // compute scattering effects
            psiNew[idx_cell] += _dE * _sigmaS[idx_energy][idx_cell] * _scatteringKernel * _psi[idx_cell];    // multiply scattering matrix with psi

            // TODO: figure out a more elegant way
            // add external source contribution
            if( _Q.size() == 1u ) {                   // constant source for all energies
                if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                    psiNew[idx_cell] += _dE * _Q[0][idx_cell][0];
                else
                    psiNew[idx_cell] += _dE * _Q[0][idx_cell];
            }
            else {
                if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                    psiNew[idx_cell] += _dE * _Q[idx_energy][idx_cell][0];
                else
                    psiNew[idx_cell] += _dE * _Q[idx_energy][idx_cell];
            }
        }
        _psi        = psiNew;
        double mass = 0.0;
        for( unsigned i = 0; i < _nCells; ++i ) {
            fluxNew[i]       = dot( _psi[i], _weights );
            _solverOutput[i] = fluxNew[i];
            mass += _areas[i] * fluxNew[i];
        }
        Save( idx_energy );
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e} {:01.5e}", _energies[idx_energy], dFlux, mass );
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
}

void SNSolver::Save( int currEnergy ) const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<std::vector<double>> scalarField( 1, _solverOutput );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), results, fieldNames, _mesh );
}
