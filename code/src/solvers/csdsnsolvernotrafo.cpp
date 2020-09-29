#include "solvers/csdsnsolvernotrafo.h"
#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/problembase.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

CSDSNSolverNoTrafo::CSDSNSolverNoTrafo( Config* settings ) : SNSolver( settings ) {
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

    _sigmaSE = std::vector<Matrix>( _energies.size(), Matrix( muMatrix.rows(), muMatrix.columns(), 0.0 ) );
    Vector angleVec( muMatrix.columns() * muMatrix.rows() );
    // store Matrix with mu values in vector format to call GetScatteringXS
    for( unsigned i = 0; i < muMatrix.rows(); ++i ) {
        for( unsigned j = 0; j < muMatrix.columns(); ++j ) {
            angleVec[i * muMatrix.columns() + j] = std::fabs( muMatrix( i, j ) );    // icru xs go from 0 to 1 due to symmetry
        }
    }

    ICRU database( angleVec, _energies );
    Matrix total;
    database.GetAngularScatteringXS( total, _sigmaTE );

    // rearrange output to matrix format
    for( unsigned n = 0; n < _energies.size(); ++n ) {
        for( unsigned i = 0; i < muMatrix.rows(); ++i ) {
            for( unsigned j = 0; j < muMatrix.columns(); ++j ) {
                _sigmaSE[n]( i, j ) = total( i * muMatrix.columns() + j, n );
            }
        }
    }

    /*
    // compute scaling s.t. scattering kernel integrates to one for chosen quadrature
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        for( unsigned p = 0; p < _nq; ++p ) {
            double cp = 0.0;
            for( unsigned q = 0; q < _nq; ++q ) {
                cp += _weights[q] * _sigmaSE[n]( p, q );
            }
            for( unsigned q = 0; q < _nq; ++q ) {
                _sigmaSE[n]( p, q ) = _sigmaTE[n] / cp * _sigmaSE[n]( p, q );
            }
        }
    }
    */

    //_s = _problem->GetStoppingPower( _energies );
    database.GetStoppingPower( _s );
    _Q = _problem->GetExternalSource( _energies );

    // recompute scattering kernel. TODO: add this to kernel function
    for( unsigned p = 0; p < _nq; ++p ) {
        for( unsigned q = 0; q < _nq; ++q ) {
            _scatteringKernel( p, q ) = 0.0;
        }
        _scatteringKernel( p, p ) = _weights[p];
    }
    // Get patient density
    _density = Vector( _nCells, 1.0 );
}

void CSDSNSolverNoTrafo::Solve() {
    std::cout << "Solve" << std::endl;
    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew( _nCells, Vector( _nq, 0.0 ) );
    double dFlux = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) log->info( "{:10}   {:10}", "E", "dFlux" );

    // revert energies and stopping power
    Vector sSave = _s;
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _energies[n] = _energies[_nEnergies - 1] - _energies[n];
        _s[n]        = sSave[_nEnergies - 1 - n];
    }

    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nq; ++k ) {
            fluxOld[j] += _weights[k] * _sol[j][k] * _s[0];
        }
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
                        psiNew[j][i] += _g->Flux( _quadPoints[i], _sol[j][i], _sol[_neighbors[j][idx_neighbor]][i], _normals[j][idx_neighbor] ) /
                                        _areas[j] / ( _density[j] * _s[n + 1] );
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[j][i] = _sol[j][i] - _dE * psiNew[j][i] - _dE * _sigmaTE[_nEnergies - n - 1] * _sol[j][i] / ( _density[j] * _s[n + 1] ) +
                               ( _s[n] / _s[n + 1] - 1 ) * _sol[j][i];
            }
            // compute scattering effects (_scatteringKernel is simply multiplication with quad weights)
            psiNew[j] += _dE * ( _sigmaSE[_nEnergies - n - 1] * _scatteringKernel * _sol[j] ) /
                         ( _density[j] * _s[n + 1] );    // multiply scattering matrix with psi
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
            fluxNew[j] = 0.0;
            for( unsigned k = 0; k < _nq; ++k ) {
                fluxNew[j] += _weights[k] * _sol[j][k] * _s[n + 1];
            }
            _dose[j] += 0.5 * _dE * ( fluxNew[j] + fluxOld[j] );    // update dose with trapezoidal rule
            _solverOutput[j] = fluxNew[j];
        }

        // Save( n );
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}  {:01.5e}  {:01.5e}", _energies[n], _dE / densityMin, dFlux );
        if( std::isinf( dFlux ) || std::isnan( dFlux ) ) break;
    }
}

void CSDSNSolverNoTrafo::Save() const {
    std::vector<std::string> fieldNames{ "dose", "normalized dose" };
    std::vector<std::vector<double>> dose( 1, _dose );
    std::vector<std::vector<double>> normalizedDose( 1, _dose );
    double maxDose = *std::max_element( _dose.begin(), _dose.end() );
    for( unsigned i = 0; i < _dose.size(); ++i ) normalizedDose[0][i] /= maxDose;
    std::vector<std::vector<std::vector<double>>> results{ dose, normalizedDose };
    ExportVTK( _settings->GetOutputFile(), results, fieldNames, _mesh );
}

void CSDSNSolverNoTrafo::Save( int currEnergy ) const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<std::vector<double>> scalarField( 1, _solverOutput );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), results, fieldNames, _mesh );
}
