#include "solvers/CSDPNSolver.h"
#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/icru.h"
#include "problems/problembase.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

CSDPNSolver::CSDPNSolver( Config* settings ) : PNSolver( settings ) {
}

void CSDPNSolver::Solve() {
}

void CSDPNSolver::PrepareVolumeOutput() {
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

void CSDPNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
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
