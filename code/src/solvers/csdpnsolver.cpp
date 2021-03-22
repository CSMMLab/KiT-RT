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
    _energies = Vector( _nEnergies, 0.0 );
    _angle    = _energies;
    _sigmaSE  = { Matrix( _nEnergies, 0.0 ) };
    _sigmaTE  = Vector( _nEnergies, 0.0 );

    Vector pos_beam = Vector{ 0.5, 0.5 };
    Vector f( _nCells );
    for( unsigned i = 0; i < _nCells; ++i ) {
        double x            = _cellMidPoints[i][0];
        double y            = _cellMidPoints[i][1];
        const double stddev = .005;
        f[i]                = normpdf( x, pos_beam[0], stddev ) * normpdf( y, pos_beam[1], stddev );
    }
    Vector IC = f * StarMAPmoments;

    Matrix sigma_t( _energies.size(), sigma_t.rows() );
    for( unsigned i = 0; i < _nSystem; ++i ) {
        Vector xs_m = blaze::column( sigma_ref, i );
        Interpolation interp( E_ref, xs_m );
        // blaze::column( sigma_t, i ) = interp( _energies );
    }

    Interpolation interpS( E_ref, S_tab );
    Vector S = interpS( _energies );
}

void CSDPNSolver::IterPreprocessing( unsigned /*idx_iter*/ ) {
    // Nothing to preprocess for PNSolver
}

void CSDPNSolver::IterPostprocessing( unsigned /*idx_iter*/ ) {
    // --- Update Solution ---
    _sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();
}

void CSDPNSolver::ComputeRadFlux() {
    double firstMomentScaleFactor = sqrt( 4 * M_PI );
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        _fluxNew[idx_cell] = _sol[idx_cell][0] * firstMomentScaleFactor;
    }
}

void CSDPNSolver::FluxUpdate() {
    // if( _reconsOrder > 1 ) {
    //    _mesh->ReconstructSlopesU( _nSystem, _solDx, _solDy, _sol );    // unstructured reconstruction
    //    //_mesh->ComputeSlopes( _nTotalEntries, _solDx, _solDy, _sol );    // unstructured reconstruction
    //}
    //// Vector solL( _nTotalEntries );
    //// Vector solR( _nTotalEntries );
    // auto solL = _sol[2];
    // auto solR = _sol[2];
    //
    //// Loop over all spatial cells
    // for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
    //
    //    // Dirichlet cells stay at IC, farfield assumption
    //    if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
    //
    //    // Reset temporary variable psiNew
    //    for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
    //        _solNew[idx_cell][idx_sys] = 0.0;
    //    }
    //
    //    // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
    //    for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); idx_neighbor++ ) {
    //
    //        // Compute flux contribution and store in psiNew to save memory
    //        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
    //            _solNew[idx_cell] += _g->Flux(
    //                _AxPlus, _AxMinus, _AyPlus, _AyMinus, _AzPlus, _AzMinus, _sol[idx_cell], _sol[idx_cell], _normals[idx_cell][idx_neighbor] );
    //        else {
    //            switch( _reconsOrder ) {
    //                // first order solver
    //                case 1:
    //                    _solNew[idx_cell] += _g->Flux( _AxPlus,
    //                                                   _AxMinus,
    //                                                   _AyPlus,
    //                                                   _AyMinus,
    //                                                   _AzPlus,
    //                                                   _AzMinus,
    //                                                   _sol[idx_cell],
    //                                                   _sol[_neighbors[idx_cell][idx_neighbor]],
    //                                                   _normals[idx_cell][idx_neighbor] );
    //                    break;
    //                // second order solver
    //                case 2:
    //                    // left status of interface
    //                    solL = _sol[idx_cell] + _solDx[idx_cell] * ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[idx_cell][0] )
    //                    +
    //                           _solDy[idx_cell] * ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[idx_cell][1] );
    //                    // right status of interface
    //                    solR = _sol[_neighbors[idx_cell][idx_neighbor]] +
    //                           _solDx[_neighbors[idx_cell][idx_neighbor]] *
    //                               ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][0] ) +
    //                           _solDy[_neighbors[idx_cell][idx_neighbor]] *
    //                               ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][1] );
    //                    // positivity checker (if not satisfied, deduce to first order)
    //                    // I manually turned it off here since Pn produces negative solutions essentially
    //                    // if( min(solL) < 0.0 || min(solR) < 0.0 ) {
    //                    //    solL = _sol[idx_cell];
    //                    //    solR = _sol[_neighbors[idx_cell][idx_neighbor]];
    //                    //}
    //
    //                    // flux evaluation
    //                    _solNew[idx_cell] +=
    //                        _g->Flux( _AxPlus, _AxMinus, _AyPlus, _AyMinus, _AzPlus, _AzMinus, solL, solR, _normals[idx_cell][idx_neighbor] );
    //                    break;
    //                // default: first order solver
    //                default:
    //                    _solNew[idx_cell] += _g->Flux( _AxPlus,
    //                                                   _AxMinus,
    //                                                   _AyPlus,
    //                                                   _AyMinus,
    //                                                   _AzPlus,
    //                                                   _AzMinus,
    //                                                   _sol[idx_cell],
    //                                                   _sol[_neighbors[idx_cell][idx_neighbor]],
    //                                                   _normals[idx_cell][idx_neighbor] );
    //            }
    //        }
    //    }
    //}
}

void CSDPNSolver::FVMUpdate( unsigned idx_energy ) {}

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
        // flux[idx_cell] = dot( _sol[idx_cell], _weights );
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
