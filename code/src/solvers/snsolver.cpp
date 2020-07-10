#include "solvers/snsolver.h"
#include "fluxes/numericalflux.h"
#include "io.h"
#include "mesh.h"
#include "settings/config.h"

#include <mpi.h>

SNSolver::SNSolver( Config* settings ) : Solver( settings ) {}

void SNSolver::Solve() {
    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _psi;

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
    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "dFlux" );

    // loop over energies (pseudo-time)
    for( unsigned n = 0; n < _nEnergies; ++n ) {

        // reconstruct slopes for higher order solver
        if( reconsOrder > 1 ) {
            _mesh->ReconstructSlopesU( _nq, psiDx, psiDy, _psi );    // unstructured reconstruction
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
                        psiNew[j][k] += _g->Flux( _quadPoints[k], _psi[j][k], _psi[j][k], _normals[j][l] );
                    else {
                        switch( reconsOrder ) {
                            // first order solver
                            case 1: psiNew[j][k] += _g->Flux( _quadPoints[k], _psi[j][k], _psi[_neighbors[j][l]][k], _normals[j][l] ); break;
                            // second order solver
                            case 2:
                                // left status of interface
                                psiL = _psi[j][k] + psiDx[j][k] * ( interfaceMidPoints[j][l][0] - cellMidPoints[j][0] ) +
                                       psiDy[j][k] * ( interfaceMidPoints[j][l][1] - cellMidPoints[j][1] );
                                // right status of interface
                                psiR = _psi[_neighbors[j][l]][k] +
                                       psiDx[_neighbors[j][l]][k] * ( interfaceMidPoints[j][l][0] - cellMidPoints[_neighbors[j][l]][0] ) +
                                       psiDy[_neighbors[j][l]][k] * ( interfaceMidPoints[j][l][1] - cellMidPoints[_neighbors[j][l]][1] );
                                // positivity check (if not satisfied, deduce to first order)
                                if( psiL < 0.0 || psiR < 0.0 ) {
                                    psiL = _psi[j][k];
                                    psiR = _psi[_neighbors[j][l]][k];
                                }
                                // flux evaluation
                                psiNew[j][k] += _g->Flux( _quadPoints[k], psiL, psiR, _normals[j][l] );
                                break;
                            // higher order solver
                            case 3: std::cout << "higher order is WIP" << std::endl; break;
                            // default: first order solver
                            default: psiNew[j][k] += _g->Flux( _quadPoints[k], _psi[j][k], _psi[_neighbors[j][l]][k], _normals[j][l] );
                        }
                    }
                }
                // time update angular flux with numerical flux and total scattering cross section
                psiNew[j][k] = _psi[j][k] - ( _dE / _areas[j] ) * psiNew[j][k] - _dE * _sigmaT[n][j] * _psi[j][k];
            }
            // compute scattering effects
            psiNew[j] += _dE * _sigmaS[n][j] * _scatteringKernel * _psi[j];    // multiply scattering matrix with psi

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
        _psi = psiNew;
        for( unsigned i = 0; i < _nCells; ++i ) {
            fluxNew[i]       = dot( _psi[i], _weights );
            _solverOutput[i] = fluxNew[i];
        }
        // Save( n );
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e}", _energies[n], dFlux );
    }
}

void SNSolver::Save() const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<double> flux( _nCells, 0.0 );
    for( unsigned i = 0; i < _nCells; ++i ) {
        flux[i] = dot( _psi[i], _weights );
    }
    std::vector<std::vector<double>> scalarField( 1, flux );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile(), results, fieldNames, _mesh );
}

void SNSolver::Save( int currEnergy ) const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<std::vector<double>> scalarField( 1, _solverOutput );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), results, fieldNames, _mesh );
}
