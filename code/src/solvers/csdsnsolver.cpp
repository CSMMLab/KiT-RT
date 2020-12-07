#include "solvers/csdsnsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/icru.h"
#include "problems/problembase.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

CSDSNSolver::CSDSNSolver( Config* settings ) : SNSolver( settings ) {
    std::cout << "Start CSDN Constructor\n";
    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

    // --- Set angle and energies
    _angle           = Vector( _settings->GetNQuadPoints(), 0.0 );    // my
    _energies        = Vector( _nEnergies, 0.0 );                     // equidistant
    double energyMin = 1e-1;
    double energyMax = 5e0;

    // --- Write equidistant energy grid
    _dE        = ComputeTimeStep( settings->GetCFL() );
    _nEnergies = unsigned( ( energyMax - energyMin ) / _dE );
    _energies.resize( _nEnergies );
    for( unsigned n = 0; n < _nEnergies; ++n ) {
        _energies[n] = energyMin + ( energyMax - energyMin ) / ( _nEnergies - 1 ) * n;
    }

    // --- Create Quadrature Points from Quad Rule
    /* Legacy Code -- unused: Will be deleted
    std::vector<double> buffer;
    for( auto q : _quadPoints )
        if( q[0] > 0.0 ) buffer.push_back( q[0] );
    std::sort( buffer.begin(), buffer.end() );
    Vector posMu( buffer.size(), buffer.data() );
    */

    // create 1D quadrature

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
            angleVec[i * muMatrix.columns() + j] = muMatrix( i, j );
        }
    }

    std::cout << "Here 3\n";

    ICRU database( angleVec, _energies, _settings );
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

    // compute scaling s.t. scattering kernel integrates to one for chosen quadrature
    /*
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
    // _s = Vector( _nEnergies, 1.0 );
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

    /*
        // test 1
        for( unsigned n = 0; n < _nEnergies; ++n ) {
            for( unsigned p = 0; p < _nq; ++p ) {
                double tmp = 0.0;
                for( unsigned q = 0; q < _nq; ++q ) {
                    tmp += _sigmaSE[n]( p, q ) * _scatteringKernel( q, q );
                }

                std::cout << tmp << " ";
            }
            std::cout << std::endl;
        }

        // exit( EXIT_FAILURE );

        // test scattering
        Vector ones( _nq, 1.0 );
        for( unsigned n = 0; n < _nEnergies; ++n ) {
            std::cout << ( _sigmaSE[_nEnergies - n - 1] * _scatteringKernel * ones )[0] << " " << _sigmaTE[_nEnergies - n - 1] << std::endl;
        }

        // exit( EXIT_FAILURE );
    */
    // Get patient density
    _density = std::vector<double>( _nCells, 1.0 );

    // Solver output
    PrepareVolumeOutput();

    std::cout << "End CSDN Constructor\n";
}

void CSDSNSolver::Solve() {
    std::cout << "Solve" << std::endl;
    auto log = spdlog::get( "event" );

    // save original energy field for boundary conditions
    auto energiesOrig = _energies;
    double energyMax  = 5e0;

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
    if( rank == 0 ) log->info( "{:10}  {:10} {:10}", "E", "  _dE / densityMin", "dFlux" );

    // --- Preprocessing ---

    // do substitution from psi to psiTildeHat (cf. Dissertation Kerstion Kuepper, Eq. 1.23)
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nq; ++k ) {
            _sol[j][k] = _sol[j][k] * _density[j] * _s[_nEnergies - 1];    // note that _s[_nEnergies - 1] is stopping power at highest energy
        }
    }

    // store transformed energies ETilde instead of E in _energies vector (cf. Dissertation Kerstion Kuepper, Eq. 1.25)
    double tmp   = 0.0;
    _energies[0] = 0.0;
    for( unsigned n = 1; n < _nEnergies; ++n ) {
        tmp          = tmp + _dE * 0.5 * ( 1.0 / _s[n] + 1.0 / _s[n - 1] );
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

    // --- end Preprocessing ---

    // loop over energies (pseudo-time)
    for( unsigned n = 0; n < _nEnergies - 1; ++n ) {
        _dE = fabs( _energies[n + 1] - _energies[n] );    // is the sign correct here?

        // --- Set Dirichlet Boundary value ---
        for( unsigned k = 0; k < _nq; ++k ) {
            if( _quadPoints[k][0] > 0 ) {
                _sol[0][k] = 1e5 * exp( -200.0 * pow( 1.0 - _quadPoints[k][0], 2 ) ) *
                             exp( -50.0 * pow( energyMax - energiesOrig[_nEnergies - n - 1], 2 ) ) * _density[0] * _s[_nEnergies - n - 1];
                if( _sol[0][k] < 0 ) {
                    ErrorMessages::Error( "boundary negative", CURRENT_FUNCTION );
                }
            }
        }

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
                psiNew[j][i] = _sol[j][i] - _dE * psiNew[j][i] - _dE * _sigmaTE[_nEnergies - n - 1] * _sol[j][i];
            }
            // compute scattering effects (_scatteringKernel is simply multiplication with quad weights)
            // std::cout << psiNew[j][0] << " | " << _sol[j][0] << std::endl;
            psiNew[j] += _dE * ( _sigmaSE[_nEnergies - n - 1] * _scatteringKernel * _sol[j] );    // multiply scattering matrix with psi
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

        for( unsigned j = 1; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            _sol[j] = psiNew[j];
        }

        // --- Postprocessing, Screen Output ---
        for( unsigned j = 0; j < _nCells; ++j ) {
            fluxNew[j] = dot( psiNew[j], _weights );
            _dose[j] += 0.5 * _dE * ( fluxNew[j] * _s[_nEnergies - n - 1] + fluxOld[j] * _s[_nEnergies - n - 2] ) /
                        _density[j];    // update dose with trapezoidal rule
            _solverOutput[j] = fluxNew[j];
        }
        // std::cout << "dose at boundary: " << _dose[0] << " \n| " << _sol[0] << std::endl;

        // --- VTK and CSV Output ---
        WriteVolumeOutput( n );
        PrintVolumeOutput( n );

        // --- Screen Output ---

        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}  {:01.5e}  {:01.5e}", _energies[n], _dE / densityMin, dFlux );
        if( std::isinf( dFlux ) || std::isnan( dFlux ) ) break;
    }
}

void CSDSNSolver::PrepareVolumeOutput() {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();

    _outputFieldNames.resize( nGroups );
    _outputFields.resize( nGroups );

    // Prepare all OutputGroups ==> Specified in option VOLUME_OUTPUT
    for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
        // Prepare all Output Fields per group

        // Different procedure, depending on the Group...
        switch( _settings->GetVolumeOutput()[idx_group] ) {
            case MINIMAL:
                // Currently only one entry ==> rad flux
                _outputFields[idx_group].resize( 1 );
                _outputFieldNames[idx_group].resize( 1 );
                _outputFields[idx_group][0].resize( _nCells );
                _outputFieldNames[idx_group][0] = "radiation flux density";
                break;

            case MEDICAL:
                // one entry per cell
                _outputFields[idx_group].resize( 1 );
                _outputFieldNames[idx_group].resize( 1 );
                _outputFields[idx_group][0].resize( _nCells );
                _outputFieldNames[idx_group][0] = std::string( "raditation dose" );
                break;

            default: ErrorMessages::Error( "Volume Output Group not defined for CSD SN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

void CSDSNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
    double mass      = 0.0;
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();

    // Compute total "mass" of the system ==> to check conservation properties
    std::vector<double> flux( _nCells, 0.0 );
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        flux[idx_cell] = dot( _sol[idx_cell], _weights );
        mass += flux[idx_cell] * _areas[idx_cell];
    }

    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_pseudoTime % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_pseudoTime == _nEnergies - 1 ) /* need sol at last iteration */ ) {

        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = flux[idx_cell];
                    }
                    break;

                case MEDICAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _dose[idx_cell];    // Remove DOSE
                    }
                    break;

                default: ErrorMessages::Error( "Volume Output Group not defined for CSD SN Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}
