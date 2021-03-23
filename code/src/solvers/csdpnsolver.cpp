#include "solvers/csdpnsolver.h"
#include "common/config.h"
#include "common/globalconstants.h"
#include "common/io.h"
#include "common/mesh.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/icru.h"
#include "problems/problembase.h"
#include "solvers/csdpn_starmap_constants.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

double normpdf( double x, double mu, double sigma ) {
    return INV_SQRT_2PI / sigma * std::exp( -( ( x - mu ) * ( x - mu ) ) / ( 2.0 * sigma * sigma ) );
}

Vector Time2Energy( const Vector& t, const double E_CutOff ) {
    Interpolation interp( E_trans, E_tab );
    Interpolation interp2( E_tab, E_trans );
    return blaze::max( 0, interp( interp2( E_CutOff, 0 ) - t ) );
}

Vector Energy2Time( const Vector& E, const double E_CutOff ) {
    Interpolation interp( E_tab, E_trans );
    return blaze::max( 0, interp( E_CutOff - E ) );
}

CSDPNSolver::CSDPNSolver( Config* settings ) : PNSolver( settings ) {
    blaze::transpose( sigma_ref );

    _dose = std::vector<double>( _settings->GetNCells(), 0.0 );

    // only dummies for compilation
    //
    //_energies = Vector( _nEnergies, 0.0 );
    //_angle    = _energies;
    //
    //_sigmaSE  = { Matrix( _nEnergies, 0.0 ) };
    //_sigmaTE  = Vector( _nEnergies, 0.0 );
    //

    Vector pos_beam = Vector{ 0.5, 0.5 };
    VectorVector IC( _nCells, Vector( _nSystem ) );
    for( unsigned i = 0; i < _nCells; ++i ) {
        double x            = _cellMidPoints[i][0];
        double y            = _cellMidPoints[i][1];
        const double stddev = .005;
        double f            = normpdf( x, pos_beam[0], stddev ) * normpdf( y, pos_beam[1], stddev );
        for( unsigned j = 0; j < _nSystem; j++ ) {
            IC[i][j] = f * StarMAPmoments[j];    // must be VectorVector
        }
    }

    Matrix sigma_t( _energies.size(), sigma_ref.rows() );
    for( unsigned i = 0; i < _nSystem; ++i ) {
        Vector xs_m = blaze::column( sigma_ref, i );
        Interpolation interp( E_ref, xs_m );
        // blaze::column( sigma_t, i ) = interp( _energies );
    }

    Interpolation interpS( E_ref, S_tab );
    _s = interpS( _energies );
}

void CSDPNSolver::SolverPreprocessing() {
    // TODO
}

void CSDPNSolver::IterPreprocessing( unsigned /*idx_iter*/ ) {
    // TODO
}

void CSDPNSolver::IterPostprocessing( unsigned idx_iter ) {
    // --- Update Solution ---
    _sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();

    unsigned n = idx_iter;
    // -- Compute Dose
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( n > 0 ) {
            _dose[j] += 0.5 * _dE * ( _fluxNew[j] * _s[_nEnergies - n - 1] + _flux[j] * _s[_nEnergies - n] ) /
                        _density[j];    // update dose with trapezoidal rule
        }
        else {
            _dose[j] += _dE * _fluxNew[j] * _s[_nEnergies - n - 1] / _density[j];
        }
        _solverOutput[j] = _fluxNew[j];    // Carefull here
        _flux[j]         = _fluxNew[j];    // Carefull here
    }
}

void CSDPNSolver::FluxUpdate() {
    // _mesh->ReconstructSlopesU( _nSystem, _solDx, _solDy, _sol );

#pragma omp parallel for
    // Loop over all spatial cells
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {

        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;

        // Reset temporary variable psiNew
        for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
            _solNew[idx_cell][idx_sys] = 0.0;
        }

        // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
        for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); idx_neighbor++ ) {

            // Compute flux contribution and store in psiNew to save memory
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                _solNew[idx_cell] += _g->Flux(
                    _AxPlus, _AxMinus, _AyPlus, _AyMinus, _AzPlus, _AzMinus, _sol[idx_cell], _sol[idx_cell], _normals[idx_cell][idx_neighbor] );
            else {

                // first order solver
                _solNew[idx_cell] += _g->Flux( _AxPlus,
                                               _AxMinus,
                                               _AyPlus,
                                               _AyMinus,
                                               _AzPlus,
                                               _AzMinus,
                                               _sol[idx_cell] / _density[idx_cell],
                                               _sol[_neighbors[idx_cell][idx_neighbor]] / _density[_neighbors[idx_cell][idx_neighbor]],
                                               _normals[idx_cell][idx_neighbor] );
            }
        }
    }
}

void CSDPNSolver::FVMUpdate( unsigned idx_energy ) {
// loop over all spatial cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // Flux update
        for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
            _solNew[idx_cell][idx_sys] = _sol[idx_cell][idx_sys] - ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_sys] /* cell averaged flux */
                                         - _dE * _sol[idx_cell][idx_sys] *
                                               ( _sigmaT[idx_energy][idx_cell]                                 /* absorbtion influence */
                                                 + _sigmaS[idx_energy][idx_cell] * _scatterMatDiag[idx_sys] ); /* scattering influence */
        }
        // Source Term
        _solNew[idx_cell][0] += _dE * _Q[0][idx_cell][0];
    }
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
                _outputFields[idx_group].resize( 2 );
                _outputFieldNames[idx_group].resize( 2 );

                // Dose
                _outputFields[idx_group][0].resize( _nCells );
                _outputFieldNames[idx_group][0] = "dose";
                // Normalized Dose
                _outputFields[idx_group][1].resize( _nCells );
                _outputFieldNames[idx_group][1] = "normalized dose";
                break;

            default: ErrorMessages::Error( "Volume Output Group not defined for CSD_SN_FP_TRAFO Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

void CSDPNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();
    double maxDose;
    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_pseudoTime % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_pseudoTime == _maxIter - 1 ) /* need sol at last iteration */ ) {

        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _fluxNew[idx_cell];
                    }
                    break;

                case MEDICAL:
                    // Compute Dose
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] += _dose[idx_cell];
                    }
                    // Compute normalized dose
                    _outputFields[idx_group][1] = _outputFields[idx_group][0];

                    maxDose = *std::max_element( _outputFields[idx_group][0].begin(), _outputFields[idx_group][0].end() );

                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][1][idx_cell] /= maxDose;
                    }
                    break;

                default: ErrorMessages::Error( "Volume Output Group not defined for CSD_SN_FP_TRAFO Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}
