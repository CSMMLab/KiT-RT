#include "solvers/csdsnsolver.h"

CSDSNSolver::CSDSNSolver( Config* settings ) : Solver( settings ) {}

void CSDSNSolver::Solve() {
    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _psi;
    double dFlux        = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) log->info( "{:10}   {:10}", "E", "dFlux" );

    // do substitution from psi to psiTildeHat (cf. Dissertation Kerstion Kuepper, Eq. 1.23)
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nq; ++k ) {
            _psi[j][k] = _psi[j][k] * _density[j] * _s[_nEnergies - 1];    // note that _s[_nEnergies - 1] is stopping power at highest energy
        }
    }

    // store transformed energies ETilde instead of E in _energies vector (cf. Dissertation Kerstion Kuepper, Eq. 1.25)
    double tmp = 0.0;
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        tmp          = tmp + _dE / _s[n];
        _energies[n] = tmp;
    }

    // store transformed energies ETildeTilde instead of ETilde in _energies vector (cf. Dissertation Kerstion Kuepper, Eq. 1.25)
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _energies[n] = _energies[_nEnergies - 1] - _energies[n];
    }

    // loop over energies (pseudo-time)
    for( unsigned n = 1; n < _nEnergies; ++n ) {
        _dE = _energies[n] - _energies[n - 1];
        // loop over all spatial cells
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            // loop over all ordinates
            for( unsigned k = 0; k < _nq; ++k ) {
                psiNew[j][k] = 0.0;
                // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                for( unsigned l = 0; l < _neighbors[j].size(); ++l ) {
                    // store flux contribution on psiNew_sigmaS to save memory
                    if( _boundaryCells[j] == BOUNDARY_TYPE::NEUMANN && _neighbors[j][l] == _nCells )
                        psiNew[j][k] += _g->Flux( _quadPoints[k] / _density[j], _psi[j][k], _psi[j][k], _normals[j][l] );
                    else
                        psiNew[j][k] += _g->Flux( _quadPoints[k] / _density[j], _psi[j][k], _psi[_neighbors[j][l]][k], _normals[j][l] );
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[j][k] = _psi[j][k] - ( _dE / _areas[j] ) * psiNew[j][k] - _dE * _sigmaT[n][j] * _psi[j][k];
            }
            // compute scattering effects
            psiNew[j] += _dE * _sigmaS[n][j] * _scatteringKernel * _psi[j];    // multiply scattering matrix with psi

            // TODO: figure out a more elegant way
            // add external source contribution
            if( _Q.size() == 1u ) {            // constant source for all energies
                if( _Q[0][j].size() == 1u )    // isotropic source
                    psiNew[j] += _dE * _Q[0][j][0] * _s[n];
                else
                    psiNew[j] += _dE * _Q[0][j] * _s[n];
            }
            else {
                if( _Q[0][j].size() == 1u )    // isotropic source
                    psiNew[j] += _dE * _Q[n][j][0] * _s[n];
                else
                    psiNew[j] += _dE * _Q[n][j] * _s[n];
            }
        }
        _psi = psiNew;

        // do backsubstitution from psiTildeHat to psi (cf. Dissertation Kerstion Kuepper, Eq. 1.23)
        for( unsigned j = 0; j < _nCells; ++j ) {
            for( unsigned k = 0; k < _nq; ++k ) {
                psiNew[j][k] = _psi[j][k] * _density[j] * _s[0];    // note that _s[0] is stopping power at highest energy
            }
        }

        for( unsigned i = 0; i < _nCells; ++i ) {
            fluxNew[i]       = dot( psiNew[i], _weights );
            _solverOutput[i] = fluxNew[i];
        }
        Save( n );
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e}", _energies[n], dFlux );
    }
}

void CSDSNSolver::Save() const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<double> flux( _nCells, 0.0 );
    for( unsigned i = 0; i < _nCells; ++i ) {
        flux[i] = dot( _psi[i], _weights );
    }
    std::vector<std::vector<double>> scalarField( 1, flux );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile(), results, fieldNames, _mesh );
}

void CSDSNSolver::Save( int currEnergy ) const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<std::vector<double>> scalarField( 1, _solverOutput );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), results, fieldNames, _mesh );
}
