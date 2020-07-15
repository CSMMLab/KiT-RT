#include "solvers/csdsnsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "physics.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

CSDSNSolver::CSDSNSolver( Config* settings ) : SNSolver( settings ) {
    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

    // Set angle and energies
    _angle    = Vector( _settings->GetNQuadPoints(), 0.0 );    // my
    _energies = Vector( _nEnergies, 0.0 );                     // equidistant
    // TODO: write meaningfull values for them!

    _sigmaS = _physics->GetScatteringXS( _energies, _angle );

    // Get patient density
    _density = Vector( _nCells, 0.0 );
}

void CSDSNSolver::Solve() {
    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _sol;
    double dFlux        = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );
    for( unsigned j = 0; j < _nCells; ++j ) {
        fluxOld[j] = dot( _sol[j], _weights );
    }
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) log->info( "{:10}   {:10}", "E", "dFlux" );

    // do substitution from psi to psiTildeHat (cf. Dissertation Kerstion Kuepper, Eq. 1.23)
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nq; ++k ) {
            _sol[j][k] = _sol[j][k] * _density[j] * _s[_nEnergies - 1];    // note that _s[_nEnergies - 1] is stopping power at highest energy
        }
    }

    // store transformed energies ETilde instead of E in _energies vector (cf. Dissertation Kerstion Kuepper, Eq. 1.25)
    double tmp = 0.0;
    for( unsigned idx_energy = 0; idx_energy < _nEnergies; ++idx_energy ) {
        tmp                   = tmp + _dE / _s[idx_energy];
        _energies[idx_energy] = tmp;
    }

    // store transformed energies ETildeTilde instead of ETilde in _energies vector (cf. Dissertation Kerstion Kuepper, Eq. 1.25)
    for( unsigned idx_energy = 0; idx_energy < _nEnergies; ++idx_energy ) {
        _energies[idx_energy] = _energies[_nEnergies - 1] - _energies[idx_energy];
    }

    // loop over energies (pseudo-time)
    for( unsigned idx_energy = 1; idx_energy < _nEnergies; ++idx_energy ) {
        _dE = fabs( _energies[idx_energy] - _energies[idx_energy - 1] );    // is the sign correct here?
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
                    // psiNew[idx_cell][idx_ord] += _g->Flux( _quadPoints[idx_ord] / _density[idx_cell],
                    //                                _psi[idx_cell][idx_ord],
                    //                              _psi[idx_cell][idx_ord],
                    //                            _normals[idx_cell][idx_neighbor] );
                    else
                        psiNew[idx_cell][idx_ord] += _g->Flux( _quadPoints[idx_ord] / _density[idx_cell],
                                                               _sol[idx_cell][idx_ord],
                                                               _sol[_neighbors[idx_cell][idx_neighbor]][idx_ord],
                                                               _normals[idx_cell][idx_neighbor] );
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[idx_cell][idx_ord] = _sol[idx_cell][idx_ord] - ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_ord] -
                                            _dE * _sigmaT[idx_energy][idx_cell] * _sol[idx_cell][idx_ord];
            }
            // compute scattering effects
            psiNew[idx_cell] += _dE * _sigmaS[idx_energy][idx_cell] * _scatteringKernel * _sol[idx_cell];    // multiply scattering matrix with psi

            // TODO: figure out a more elegant way
            // add external source contribution
            if( _Q.size() == 1u ) {                   // constant source for all energies
                if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                    psiNew[idx_cell] += _dE * _Q[0][idx_cell][0] * _s[_nEnergies - idx_energy - 1];
                else
                    psiNew[idx_cell] += _dE * _Q[0][idx_cell] * _s[_nEnergies - idx_energy - 1];
            }
            else {
                if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                    psiNew[idx_cell] += _dE * _Q[idx_energy][idx_cell][0] * _s[_nEnergies - idx_energy - 1];
                else
                    psiNew[idx_cell] += _dE * _Q[idx_energy][idx_cell] * _s[_nEnergies - idx_energy - 1];
            }
        }
        _sol = psiNew;

        // do backsubstitution from psiTildeHat to psi (cf. Dissertation Kerstion Kuepper, Eq. 1.23)
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
            for( unsigned idx_ord = 0; idx_ord < _nq; ++idx_ord ) {
                psiNew[idx_cell][idx_ord] = _sol[idx_cell][idx_ord] * _density[idx_cell] *
                                            _s[_nEnergies - idx_energy - 1];    // note that _s[0] is stopping power at lowest energy
            }
        }

        for( unsigned j = 0; j < _nCells; ++j ) {
            fluxNew[j] = dot( psiNew[j], _weights );
            _dose[j] += 0.5 * _dE * ( fluxNew[j] + fluxOld[j] ) / _density[j];    // update dose with trapezoidal rule
            _solverOutput[j] = fluxNew[j];
        }

        Save( idx_energy );
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e}", _energies[idx_energy], dFlux );
    }
}

void CSDSNSolver::Save() const {
    std::vector<std::string> fieldNames{ "dose" };
    std::vector<std::vector<double>> scalarField( 1, _dose );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile(), results, fieldNames, _mesh );
}

void CSDSNSolver::Save( int currEnergy ) const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<std::vector<double>> scalarField( 1, _solverOutput );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), results, fieldNames, _mesh );
}
