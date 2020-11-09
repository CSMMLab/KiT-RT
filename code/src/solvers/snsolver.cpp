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
    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _sol;

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
    for( unsigned i = 0; i < _nCells; ++i ) {
        for( unsigned k = 0; k < _mesh->GetDim(); ++k ) {
            for( unsigned j = 0; j < _neighbors[i].size() - 1; ++j ) {
                interfaceMidPoints[i][j][k] = 0.5 * ( nodes[cells[i][j]][k] + nodes[cells[i][j + 1]][k] );
            }
            interfaceMidPoints[i][_neighbors[i].size() - 1][k] = 0.5 * ( nodes[cells[i][_neighbors[i].size() - 1]][k] + nodes[cells[i][0]][k] );
        }
    }

    // distance between cell center to interface center
    VectorVector cellDx( _nCells, Vector( _mesh->GetNumNodesPerCell(), 1e-10 ) );
    VectorVector cellDy( _nCells, Vector( _mesh->GetNumNodesPerCell(), 1e-10 ) );
    for( unsigned i = 0; i < _nCells; ++i ) {
        for( unsigned j = 0; j < _mesh->GetNumNodesPerCell(); ++j ) {
            cellDx[i][j] = interfaceMidPoints[i][j][0] - cellMidPoints[i][0];
            cellDy[i][j] = interfaceMidPoints[i][j][1] - cellMidPoints[j][1];
        }
    }

    double dFlux = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "dFlux", "mass" );

    // if( rank == 0 ) log->info( "{:01.5e}   {:01.5e} {:01.5e}", 0.0, dFlux, mass1 );

    // loop over energies (pseudo-time)
    for( unsigned n = 0; n < _nEnergies; ++n ) {

        // reconstruct slopes for higher order solver
        if( reconsOrder > 1 ) {
            _mesh->ReconstructSlopesU( _nq, psiDx, psiDy, _sol );    // unstructured reconstruction
            //_mesh->ReconstructSlopesS( _nq, psiDx, psiDy, _psi );    // structured reconstruction (not stable currently)
        }

        // loop over all spatial cells
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) continue;
            // loop over all ordinates
            for( unsigned k = 0; k < _nq; ++k ) {
                psiNew[j][k] = 0.0;
                // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                for( unsigned l = 0; l < _neighbors[j].size(); ++l ) {
                    // store flux contribution on psiNew_sigmaS to save memory
                    if( _boundaryCells[j] == BOUNDARY_TYPE::NEUMANN && _neighbors[j][l] == _nCells )
                        psiNew[j][k] += _g->Flux( _quadPoints[k], _sol[j][k], _sol[j][k], _normals[j][l] );
                    else {
                        switch( reconsOrder ) {
                            // first order solver
                            case 1: psiNew[j][k] += _g->Flux( _quadPoints[k], _sol[j][k], _sol[_neighbors[j][l]][k], _normals[j][l] ); break;
                            // second order solver
                            case 2:
                                // left status of interface
                                psiL = _sol[j][k] + psiDx[j][k] * ( interfaceMidPoints[j][l][0] - cellMidPoints[j][0] ) +
                                       psiDy[j][k] * ( interfaceMidPoints[j][l][1] - cellMidPoints[j][1] );
                                // right status of interface
                                psiR = _sol[_neighbors[j][l]][k] +
                                       psiDx[_neighbors[j][l]][k] * ( interfaceMidPoints[j][l][0] - cellMidPoints[_neighbors[j][l]][0] ) +
                                       psiDy[_neighbors[j][l]][k] * ( interfaceMidPoints[j][l][1] - cellMidPoints[_neighbors[j][l]][1] );
                                // positivity check (if not satisfied, deduce to first order)
                                if( psiL < 0.0 || psiR < 0.0 ) {
                                    psiL = _sol[j][k];
                                    psiR = _sol[_neighbors[j][l]][k];
                                }
                                // flux evaluation
                                psiNew[j][k] += _g->Flux( _quadPoints[k], psiL, psiR, _normals[j][l] );
                                break;
                            // higher order solver
                            case 3: std::cout << "higher order is WIP" << std::endl; break;
                            // default: first order solver
                            default: psiNew[j][k] += _g->Flux( _quadPoints[k], _sol[j][k], _sol[_neighbors[j][l]][k], _normals[j][l] );
                        }
                    }
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[j][k] = _sol[j][k] - ( _dE / _areas[j] ) * psiNew[j][k] - _dE * _sigmaT[n][j] * _sol[j][k];
            }
            // compute scattering effects
            psiNew[j] += _dE * _sigmaS[n][j] * _scatteringKernel * _sol[j];    // multiply scattering matrix with psi

            // TODO: figure out a more elegant way
            // add external source contribution
            if( _Q.size() == 1u ) {            // constant source for all energies
                if( _Q[0][j].size() == 1u )    // isotropic source
                    psiNew[j] += _dE * _Q[0][j][0];
                else
                    psiNew[j] += _dE * _Q[0][j];
            }
            else {
                if( _Q[0][j].size() == 1u )    // isotropic source
                    psiNew[j] += _dE * _Q[n][j][0];
                else
                    psiNew[j] += _dE * _Q[n][j];
            }
        }
        _sol = psiNew;
        for( unsigned i = 0; i < _nCells; ++i ) {
            fluxNew[i] = dot( _sol[i], _weights );
        }

        // --- VTK and CSV Output ---
        mass = WriteOutputFields( n );
        Save( n );

        // --- Screen Output ---
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e} {:01.5e}", _energies[n], dFlux, mass );
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
