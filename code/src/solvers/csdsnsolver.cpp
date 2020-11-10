#include "solvers/csdsnsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/problembase.h"
#include "toolboxes/errormessages.h"

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
    for( unsigned k = 0; k < _settings->GetNQuadPoints(); ++k ) {
        _angle[k] = _quadPoints[k][2];
    }

    _sigmaSE = _problem->GetScatteringXSE( _energies, _angle );
    _sigmaTE = _problem->GetTotalXSE( _energies );
    _s       = _problem->GetStoppingPower( _energies );
    _Q       = _problem->GetExternalSource( _energies );

    // Get patient density
    _density = Vector( _nCells, 1.0 );

    // Solver output
    PrepareOutputFields();
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
                                                  _normals[j][idx_neighbor] );
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[j][i] = _sol[j][i] - ( _dE / _areas[j] ) * psiNew[j][i] - _dE * _sigmaTE[n] * _sol[j][i];
            }
            // compute scattering effects (_scatteringKernel is simply multiplication with quad weights)
            psiNew[j] += _dE * ( blaze::trans( _sigmaSE[n] ) * _scatteringKernel * _sol[j] );    // multiply scattering matrix with psi
            // TODO: Check if _sigmaS^T*psi is correct

            // TODO: figure out a more elegant way
            // add external source contribution
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
        }
        _sol = psiNew;
        // do backsubstitution from psiTildeHat to psi (cf. Dissertation Kerstion Kuepper, Eq. 1.23)
        for( unsigned j = 0; j < _nCells; ++j ) {
            for( unsigned i = 0; i < _nq; ++i ) {
                psiNew[j][i] = _sol[j][i] * _density[j] * _s[_nEnergies - n - 1];    // note that _s[0] is stopping power at lowest energy
            }
        }

        for( unsigned j = 0; j < _nCells; ++j ) {
            fluxNew[j] = dot( psiNew[j], _weights );
            _dose[j] += 0.5 * _dE * ( fluxNew[j] + fluxOld[j] ) / _density[j];    // update dose with trapezoidal rule
            _solverOutput[j] = fluxNew[j];
        }

        // --- VTK and CSV Output ---
        WriteOutputFields( n );
        Save( n );

        // --- Screen Output ---
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}  {:01.5e}  {:01.5e}", _energies[n], _dE / densityMin, dFlux );
        if( std::isinf( dFlux ) || std::isnan( dFlux ) ) break;
    }
}

void CSDSNSolver::PrepareOutputFields() {
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

            case DOSE:
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

void CSDSNSolver::WriteOutputFields( unsigned idx_pseudoTime ) {
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

                case DOSE:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _dose[idx_cell];    // Remove DOSE
                    }
                    break;

                default: ErrorMessages::Error( "Volume Output Group not defined for CSD SN Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}
