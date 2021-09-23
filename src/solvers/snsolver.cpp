#include "solvers/snsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"
#include "toolboxes/errormessages.h"

// externals
#include "spdlog/spdlog.h"
#include <iostream>
#include <mpi.h>

SNSolver::SNSolver( Config* settings ) : SolverBase( settings ) {
    _psiDx = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    _psiDy = VectorVector( _nCells, Vector( _nq, 0.0 ) );

    _quadPoints = _quadrature->GetPoints();
    _weights    = _quadrature->GetWeights();

    ScatteringKernel* k = ScatteringKernel::CreateScatteringKernel( settings->GetKernelName(), _quadrature );
    _scatteringKernel   = k->GetScatteringKernel();
    delete k;
}

void SNSolver::IterPreprocessing( unsigned /*idx_iter*/ ) {
    // Nothing to do for SNSolver
}

void SNSolver::IterPostprocessing( unsigned /*idx_iter*/ ) {
    // --- Update Solution ---
    _sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();
}

void SNSolver::ComputeRadFlux() {
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        _fluxNew[idx_cell] = blaze::dot( _sol[idx_cell], _weights );
    }
}

void SNSolver::FluxUpdate() {

    // Legacy Code form second order reconstruction:

    /*
    // reconstruction order
    unsigned reconsOrder = _settings->GetReconsOrder();

    // left and right angular flux of interface, used in numerical flux evaluation
    double psiL;
    double psiR;

    // derivatives of angular flux in x and y directions
    VectorVector psiDx( _nCells, Vector( _nq, 0.0 ) );
    VectorVector psiDy( _nCells, Vector( _nq, 0.0 ) );


    // geometric variables for derivatives computation
    auto nodes         = _mesh->GetNodes();
    auto cells         = _mesh->GetCells();
    auto cellMidPoints = _mesh->GetCellMidPoints();


    // center location of cell interfaces
    std::vector<std::vector<Vector>> interfaceMidPoints( _nCells, std::vector<Vector>( _mesh->GetNumNodesPerCell(), Vector( 2, 1e-10 ) ) );
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        for( unsigned k = 0; k < _mesh->GetDim(); ++k ) {
            for( unsigned j = 0; j < _neighbors[idx_cell].size() - 1; ++j ) {
                interfaceMidPoints[idx_cell][j][k] = 0.5 * ( nodes[cells[idx_cell][j]][k] + nodes[cells[idx_cell][j + 1]][k] );
            }
            interfaceMidPoints[idx_cell][_neighbors[idx_cell].size() - 1][k] =
                0.5 * ( nodes[cells[idx_cell][_neighbors[idx_cell].size() - 1]][k] + nodes[cells[idx_cell][0]][k] );
        }
    }

     distance between cell center to interface center
     VectorVector cellDx( _nCells, Vector( _mesh->GetNumNodesPerCell(), 1e-10 ) );
     VectorVector cellDy( _nCells, Vector( _mesh->GetNumNodesPerCell(), 1e-10 ) );

     for( unsigned i = 0; i < _nCells; ++i ) {
        for( unsigned j = 0; j < _mesh->GetNumNodesPerCell(); ++j ) {
            cellDx[i][j] = interfaceMidPoints[i][j][0] - cellMidPoints[i][0];
            cellDy[i][j] = interfaceMidPoints[i][j][1] - cellMidPoints[j][1];
        }
     }
    */
    /*
    // reconstruct slopes for higher order solver
     if( reconsOrder > 1 ) {
        _mesh->ReconstructSlopesU( _nq, psiDx, psiDy, _sol );    // unstructured reconstruction
        //_mesh->ReconstructSlopesS( _nq, psiDx, psiDy, _psi );    // structured reconstruction (not stable currently)
    }
    */
    if( _reconsOrder > 1 ) {
        //_mesh->ReconstructSlopesU( _nq, _psiDx, _psiDy, _sol );    // unstructured reconstruction
        _mesh->ReconstructSlopesS( _nq, _psiDx, _psiDy, _sol );    // structured reconstruction (not stable currently)
    }
    double psiL;
    double psiR;
    // Loop over all spatial cells
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // Loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nq; ++idx_quad ) {
            // Reset temporary variable
            _solNew[idx_cell][idx_quad] = 0.0;
            // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); ++idx_neighbor ) {
                // store flux contribution on psiNew_sigmaS to save memory
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                    _solNew[idx_cell][idx_quad] +=
                        _g->Flux( _quadPoints[idx_quad], _sol[idx_cell][idx_quad], _sol[idx_cell][idx_quad], _normals[idx_cell][idx_neighbor] );
                else {
                    switch( _reconsOrder ) {
                        // first order solver
                        case 1:
                            _solNew[idx_cell][idx_quad] += _g->Flux( _quadPoints[idx_quad],
                                                                     _sol[idx_cell][idx_quad],
                                                                     _sol[_neighbors[idx_cell][idx_neighbor]][idx_quad],
                                                                     _normals[idx_cell][idx_neighbor] );
                            break;
                        // second order solver
                        case 2:
                            // left status of interface
                            psiL = _sol[idx_cell][idx_quad] +
                                   _psiDx[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[idx_cell][0] ) +
                                   _psiDy[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[idx_cell][1] );
                            // right status of interface
                            psiR = _sol[_neighbors[idx_cell][idx_neighbor]][idx_quad] +
                                   _psiDx[_neighbors[idx_cell][idx_neighbor]][idx_quad] *
                                       ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][0] ) +
                                   _psiDy[_neighbors[idx_cell][idx_neighbor]][idx_quad] *
                                       ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][1] );
                            // positivity check (if not satisfied, deduce to first order)
                            if( psiL < 0.0 || psiR < 0.0 ) {
                                psiL = _sol[idx_cell][idx_quad];
                                psiR = _sol[_neighbors[idx_cell][idx_neighbor]][idx_quad];
                            }
                            // flux evaluation
                            _solNew[idx_cell][idx_quad] += _g->Flux( _quadPoints[idx_quad], psiL, psiR, _normals[idx_cell][idx_neighbor] );
                            break;
                        // higher order solver
                        case 3:
                            std::cout << "higher order is WIP" << std::endl;
                            break;
                            // default: first order solver
                        default:
                            _solNew[idx_cell][idx_quad] += _g->Flux( _quadPoints[idx_quad],
                                                                     _sol[idx_cell][idx_quad],
                                                                     _sol[_neighbors[idx_cell][idx_neighbor]][idx_quad],
                                                                     _normals[idx_cell][idx_neighbor] );
                    }
                }
            }
        }
    }
}

void SNSolver::FVMUpdate( unsigned idx_energy ) {
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nq; ++idx_quad ) {
            // time update angular flux with numerical flux and total scattering cross section
            _solNew[idx_cell][idx_quad] = _sol[idx_cell][idx_quad] - ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_quad] -
                                          _dE * _sigmaT[idx_energy][idx_cell] * _sol[idx_cell][idx_quad];
        }
        // compute scattering effects
        _solNew[idx_cell] += _dE * _sigmaS[idx_energy][idx_cell] * _scatteringKernel * _sol[idx_cell];    // multiply scattering matrix with psi

        // Source Term
        if( _Q.size() == 1u ) {                   // constant source for all energies
            if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                _solNew[idx_cell] += _dE * _Q[0][idx_cell][0];
            else
                _solNew[idx_cell] += _dE * _Q[0][idx_cell];
        }
        else {
            if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                _solNew[idx_cell] += _dE * _Q[idx_energy][idx_cell][0];
            else
                _solNew[idx_cell] += _dE * _Q[idx_energy][idx_cell];
        }
    }
}

void SNSolver::PrepareVolumeOutput() {
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

            case ANALYTIC:
                // one entry per cell
                _outputFields[idx_group].resize( 1 );
                _outputFieldNames[idx_group].resize( 1 );
                _outputFields[idx_group][0].resize( _nCells );
                _outputFieldNames[idx_group][0] = std::string( "analytic radiation flux density" );
                break;

            default: ErrorMessages::Error( "Volume Output Group not defined for SN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

void SNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();

    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_pseudoTime % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_pseudoTime == _nEnergies - 1 ) /* need sol at last iteration */ ) {

        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _fluxNew[idx_cell];
                    }
                    break;

                case ANALYTIC:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        double time                           = idx_pseudoTime * _dE;
                        _outputFields[idx_group][0][idx_cell] = _problem->GetAnalyticalSolution(
                            _mesh->GetCellMidPoints()[idx_cell][0], _mesh->GetCellMidPoints()[idx_cell][1], time, _sigmaS[idx_pseudoTime][idx_cell] );
                    }
                    break;

                default: ErrorMessages::Error( "Volume Output Group not defined for SN Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}
