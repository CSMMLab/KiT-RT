#include "solvers/snsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

SNSolver::SNSolver( Config* settings ) : Solver( settings ) {

    _quadPoints = _quadrature->GetPoints();
    _weights    = _quadrature->GetWeights();

    ScatteringKernel* k = ScatteringKernel::CreateScatteringKernel( settings->GetKernelName(), _quadrature );
    _scatteringKernel   = k->GetScatteringKernel();
    delete k;

    // Solver output
    PrepareOutputFields();
}

void SNSolver::Solve() {

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    auto log = spdlog::get( "event" );

    VectorVector psiNew = _sol;

    double mass  = 0.0;
    double dFlux = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );

    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "dFlux", "mass" );

    // loop over energies (pseudo-time)
    for( unsigned idx_energy = 0; idx_energy < _nEnergies; ++idx_energy ) {

        // --- Prepare Boundaries and temp variables
        IterPreprocessing();

        // --- Compute Fluxes ---
        FluxUpdate( psiNew );

        // --- Compute Finite Volume Step ---
        FVMUpdate( psiNew, idx_energy );

        // --- Update Solution ---
        _sol = psiNew;

        // --- VTK and CSV Output ---
        mass = WriteOutputFields( idx_energy );
        Save( idx_energy );

        // --- Screen Output ---
        for( unsigned i = 0; i < _nCells; ++i ) {
            fluxNew[i] = dot( _sol[i], _weights );
        }

        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e} {:01.5e}", _energies[idx_energy], dFlux, mass );
    }
}

void SNSolver::IterPreprocessing() {
    // Nothing to do for SNSolver
}

void SNSolver::FluxUpdate( VectorVector& psiNew ) {

    // Legacy Code form second order reconstruction:

    /*
    // reconstruction order
    unsigned reconsOrder = _settings->GetReconsOrder();

    // left and right angular flux of interface, used in numerical flux evaluation
    double psiL;
    double psiR;
    double mass;

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

    // Loop over all spatial cells
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // Loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nq; ++idx_quad ) {
            // Reset temporary variable
            psiNew[idx_cell][idx_quad] = 0.0;
            // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); ++idx_neighbor ) {
                // store flux contribution on psiNew_sigmaS to save memory
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                    psiNew[idx_cell][idx_quad] +=
                        _g->Flux( _quadPoints[idx_quad], _sol[idx_cell][idx_quad], _sol[idx_cell][idx_quad], _normals[idx_cell][idx_neighbor] );
                else {
                    /*                    switch( reconsOrder ) {
                                            // first order solver
                                            case 1:
                                                psiNew[idx_cells][idx_quad] += _g->Flux( _quadPoints[idx_quad],
                                                                                         _sol[idx_cells][idx_quad],
                                                                                         _sol[_neighbors[idx_cells][idx_neighbor]][idx_quad],
                                                                                         _normals[idx_cells][idx_neighbor] );
                                                break;
                                            // second order solver
                                            case 2:
                                                // left status of interface
                                                psiL = _sol[idx_cells][idx_quad] +
                                                       psiDx[idx_cells][idx_quad] * ( interfaceMidPoints[idx_cells][idx_neighbor][0] -
                       cellMidPoints[idx_cells][0] ) + psiDy[idx_cells][idx_quad] * ( interfaceMidPoints[idx_cells][idx_neighbor][1] -
                       cellMidPoints[idx_cells][1] );
                                                // right status of interface
                                                psiR = _sol[_neighbors[idx_cells][idx_neighbor]][idx_quad] +
                                                       psiDx[_neighbors[idx_cells][idx_neighbor]][idx_quad] *
                                                           ( interfaceMidPoints[idx_cells][idx_neighbor][0] -
                       cellMidPoints[_neighbors[idx_cells][idx_neighbor]][0] ) + psiDy[_neighbors[idx_cells][idx_neighbor]][idx_quad] * (
                       interfaceMidPoints[idx_cells][idx_neighbor][1] - cellMidPoints[_neighbors[idx_cells][idx_neighbor]][1] );
                                                // positivity check (if not satisfied, deduce to first order)
                                                if( psiL < 0.0 || psiR < 0.0 ) {
                                                    psiL = _sol[idx_cells][idx_quad];
                                                    psiR = _sol[_neighbors[idx_cells][idx_neighbor]][idx_quad];
                                                }
                                                // flux evaluation
                                                psiNew[idx_cells][idx_quad] += _g->Flux( _quadPoints[idx_quad], psiL, psiR,
                       _normals[idx_cells][idx_neighbor] ); break;
                                            // higher order solver
                                            case 3: std::cout << "higher order is WIP" << std::endl; break;
                                            // default: first order solver
                                            default:
                                            */
                    psiNew[idx_cell][idx_quad] += _g->Flux( _quadPoints[idx_quad],
                                                            _sol[idx_cell][idx_quad],
                                                            _sol[_neighbors[idx_cell][idx_neighbor]][idx_quad],
                                                            _normals[idx_cell][idx_neighbor] );
                    // }
                }
            }
        }
    }
}

void SNSolver::FVMUpdate( VectorVector& psiNew, unsigned idx_energy ) {
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nq; ++idx_quad ) {
            // time update angular flux with numerical flux and total scattering cross section
            psiNew[idx_cell][idx_quad] = _sol[idx_cell][idx_quad] - ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_quad] -
                                         _dE * _sigmaT[idx_energy][idx_cell] * _sol[idx_cell][idx_quad];
        }
        // compute scattering effects
        psiNew[idx_cell] += _dE * _sigmaS[idx_energy][idx_cell] * _scatteringKernel * _sol[idx_cell];    // multiply scattering matrix with psi

        // Source Term
        if( _Q.size() == 1u ) {                   // constant source for all energies
            if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                psiNew[idx_cell] += _dE * _Q[0][idx_cell][0];
            else
                psiNew[idx_cell] += _dE * _Q[0][idx_cell];
        }
        else {
            if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                psiNew[idx_cell] += _dE * _Q[idx_energy][idx_cell][0];
            else
                psiNew[idx_cell] += _dE * _Q[idx_energy][idx_cell];
        }
    }
}

void SNSolver::PrepareOutputFields() {
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

double SNSolver::WriteOutputFields( unsigned idx_pseudoTime ) {
    double mass      = 0.0;
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();

    // Compute total "mass" of the system ==> to check conservation properties
    std::vector<double> flux( _nCells, 0.0 );
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        flux[idx_cell] = dot( _sol[idx_cell], _weights );
        mass += flux[idx_cell] * _areas[idx_cell];
    }

    if( ( _settings->GetOutputFrequency() != 0 && idx_pseudoTime % (unsigned)_settings->GetOutputFrequency() == 0 ) ||
        ( idx_pseudoTime == _nEnergies - 1 ) /* need sol at last iteration */ ) {

        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = flux[idx_cell];
                    }
                    break;

                case ANALYTIC:
                    // Compute total "mass" of the system ==> to check conservation properties
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {

                        double time = idx_pseudoTime * _dE;

                        _outputFields[idx_group][0][idx_cell] = _problem->GetAnalyticalSolution(
                            _mesh->GetCellMidPoints()[idx_cell][0], _mesh->GetCellMidPoints()[idx_cell][1], time, _sigmaS[idx_pseudoTime][idx_cell] );
                    }
                    break;

                default: ErrorMessages::Error( "Volume Output Group not defined for SN Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
    return mass;
}
