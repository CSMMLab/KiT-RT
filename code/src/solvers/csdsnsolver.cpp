#include "solvers/csdsnsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/problembase.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

CSDSNSolver::CSDSNSolver( Config* settings ) : SNSolver( settings ) {
    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

    // Set angle and energies
    _angle           = Vector( _settings->GetNQuadPoints(), 0.0 );    // my
    _energies        = Vector( _nEnergies, 0.0 );                     // equidistant
    double energyMin = 1e-1;
    double energyMax = 5e0;
    // write equidistant energy grid

    _dE        = ComputeTimeStep( settings->GetCFL() );
    _nEnergies = unsigned( ( energyMax - energyMin ) / _dE );
    _energies.resize( _nEnergies );
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _energies[n] = energyMin + ( energyMax - energyMin ) / ( _nEnergies - 1 ) * n;
    }
    // write mu grid
    Matrix muMatrix( _settings->GetNQuadPoints(), _settings->GetNQuadPoints() );
    for( unsigned l = 0; l < _settings->GetNQuadPoints(); ++l ) {
        for( unsigned k = 0; k < _settings->GetNQuadPoints(); ++k ) {
            double inner = 0.0;
            for( unsigned j = 0; j < 3; ++j ) {
                inner += _quadPoints[l][j] * _quadPoints[k][j];    // compute mu at Omega_l*Omega_k
            }
            muMatrix( l, k ) = inner;
        }
    }

    _sigmaSE = _problem->GetScatteringXSE( _energies, muMatrix );
    _sigmaTE = _problem->GetTotalXSE( _energies );
    _s       = _problem->GetStoppingPower( _energies );
    _Q       = _problem->GetExternalSource( _energies );

    double dMu   = 0.00001;
    unsigned nMu = unsigned( 2 / dMu );
    Vector muGrid( nMu );
    for( unsigned n = 0; n < nMu; ++n ) {
        muGrid[n] = -1.0 + 2.0 / ( nMu - 1 ) * n;
    }
    VectorVector tmp = _problem->GetScatteringXSE( _energies, muGrid );
    double sigmaT    = 0.0;
    for( unsigned n = 0; n < nMu; ++n ) {
        sigmaT += tmp[0][n] * dMu;
    }
    std::cout << "int(sigmaS) at energy " << 2.0 * M_PI * sigmaT << std::endl;
    std::cout << "sigmaT at energy " << _sigmaTE[0] << std::endl;

    // Get patient density
    _density = Vector( _nCells, 1.0 );
}

void CSDSNSolver::Solve() {
    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew( _nCells, Vector( _nq, 0.0 ) );
    double dFlux = 1e10;
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
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        tmp          = tmp + _dE / _s[n];
        _energies[n] = tmp;
    }

    // store transformed energies ETildeTilde instead of ETilde in _energies vector (cf. Dissertation Kerstion Kuepper, Eq. 1.25)
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _energies[n] = _energies[_nEnergies - 1] - _energies[n];
    }

    // determine minimal density for CFL computation
    double densityMin = _density[0];
    for( unsigned j = 1; j < _nCells; ++j ) {
        if( densityMin > _density[j] ) densityMin = _density[j];
    }

    // cross sections do not need to be transformed to ETilde energy grid since e.g. TildeSigmaT(ETilde) = SigmaT(E(ETilde))

    // loop over energies (pseudo-time)
    for( unsigned n = 0; n < _nEnergies - 1; ++n ) {
        _dE = fabs( _energies[n + 1] - _energies[n] );    // is the sign correct here?
        // loop over all spatial cells
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            // loop over all ordinates
            for( unsigned i = 0; i < _nq; ++i ) {
                psiNew[j][i] = 0.0;
                // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[j].size(); ++idx_neighbor ) {
                    // store flux contribution on psiNew_sigmaS to save memory
                    if( _boundaryCells[j] == BOUNDARY_TYPE::NEUMANN && _neighbors[j][idx_neighbor] == _nCells )
                        continue;    // adiabatic wall, add nothing
                    else
                        psiNew[j][i] += _g->Flux( _quadPoints[i],
                                                  _sol[j][i] / _density[j],
                                                  _sol[_neighbors[j][idx_neighbor]][i] / _density[_neighbors[j][idx_neighbor]],
                                                  _normals[j][idx_neighbor] ) /
                                        _areas[j];
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[j][i] = _sol[j][i] - _dE * psiNew[j][i] - _dE * _sigmaTE[n] * _sol[j][i];
            }
            // compute scattering effects (_scatteringKernel is simply multiplication with quad weights)
            psiNew[j] += _dE * ( _sigmaSE[n] * _scatteringKernel * _sol[j] );    // multiply scattering matrix with psi
            // TODO: Check if _sigmaS^T*psi is correct

            // TODO: figure out a more elegant way
            // add external source contribution
            /*
            if( _Q.size() == 1u ) {            // constant source for all energies
                if( _Q[0][j].size() == 1u )    // isotropic source
                    psiNew[j] += _dE * _Q[0][j][0] * _s[_nEnergies - n - 1];
                else
                    psiNew[j] += _dE * _Q[0][j] * _s[_nEnergies - n - 1];
            }
            else {
                if( _Q[0][j].size() == 1u )    // isotropic source
                    psiNew[j] += _dE * _Q[n][j][0] * _s[_nEnergies - n - 1];
                else
                    psiNew[j] += _dE * _Q[n][j] * _s[_nEnergies - n - 1];
            }
            */
        }
        _sol = psiNew;

        for( unsigned j = 0; j < _nCells; ++j ) {
            fluxNew[j] = dot( psiNew[j], _weights );
            _dose[j] += 0.5 * _dE * ( fluxNew[j] / _s[_energies.size() - n - 1] + fluxOld[j] / _s[_energies.size() - n - 2] ) /
                        _density[j];    // update dose with trapezoidal rule
            _solverOutput[j] = fluxNew[j];
        }

        // Save( n );
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}  {:01.5e}  {:01.5e}", _energies[n], _dE / densityMin, dFlux );
        if( std::isinf( dFlux ) || std::isnan( dFlux ) ) break;
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
