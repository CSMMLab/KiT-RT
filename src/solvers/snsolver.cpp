#include "solvers/snsolver.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "fluxes/numericalflux.hpp"
#include "kernels/scatteringkernelbase.hpp"
#include "problems/problembase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

// externals
#include "spdlog/spdlog.h"
#include <iostream>
#include <mpi.h>

SNSolver::SNSolver( Config* settings ) : SolverBase( settings ) {
    _quadPoints = _quadrature->GetPoints();
    _weights    = _quadrature->GetWeights();

    ScatteringKernel* k = ScatteringKernel::CreateScatteringKernel( settings->GetKernelName(), _quadrature );
    _scatteringKernel   = k->GetScatteringKernel();

    // Limiter variables
    _solDx   = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    _solDy   = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    _limiter = VectorVector( _nCells, Vector( _nq, 0.0 ) );

    delete k;
}

void SNSolver::IterPreprocessing( unsigned /*idx_iter*/ ) {
    // Slope Limiter computation
    if( _reconsOrder > 1 ) {
        _mesh->ComputeSlopes( _nq, _solDx, _solDy, _sol );
        _mesh->ComputeLimiter( _nq, _solDx, _solDy, _sol, _limiter );
    }
}

void SNSolver::IterPostprocessing( unsigned /*idx_iter*/ ) {
    // --- Update Solution ---
    _sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();
}

void SNSolver::ComputeRadFlux() {
    double firstMomentScaleFactor = 4 * M_PI;
    if( _settings->GetProblemName() == PROBLEM_Aircavity1D || _settings->GetProblemName() == PROBLEM_Linesource1D ||
        _settings->GetProblemName() == PROBLEM_Checkerboard1D || _settings->GetProblemName() == PROBLEM_Meltingcube1D ) {
        firstMomentScaleFactor = 2.0;
    }
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        _fluxNew[idx_cell] = blaze::dot( _sol[idx_cell], _weights ) / firstMomentScaleFactor;
    }
}

void SNSolver::FluxUpdate() {
    if( _settings->GetProblemName() == PROBLEM_Aircavity1D || _settings->GetProblemName() == PROBLEM_Linesource1D ||
        _settings->GetProblemName() == PROBLEM_Checkerboard1D || _settings->GetProblemName() == PROBLEM_Meltingcube1D ) {
        FluxUpdatePseudo1D();
    }
    else {
        FluxUpdatePseudo2D();
    }
}

void SNSolver::FluxUpdatePseudo1D() {
// Loop over all spatial cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        double solL;
        double solR;
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // Loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nq; ++idx_quad ) {
            // Reset temporary variable
            _solNew[idx_cell][idx_quad] = 0.0;
            // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_nbr = 0; idx_nbr < _neighbors[idx_cell].size(); ++idx_nbr ) {
                // store flux contribution on psiNew_sigmaS to save memory
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_nbr] == _nCells )
                    _solNew[idx_cell][idx_quad] +=
                        _g->Flux1D( _quadPoints[idx_quad], _sol[idx_cell][idx_quad], _sol[idx_cell][idx_quad], _normals[idx_cell][idx_nbr] );
                else {
                    unsigned int nbr_glob = _neighbors[idx_cell][idx_nbr];    // global idx of neighbor cell

                    switch( _reconsOrder ) {
                        // first order solver
                        case 1:
                            _solNew[idx_cell][idx_quad] +=
                                _g->Flux1D( _quadPoints[idx_quad], _sol[idx_cell][idx_quad], _sol[nbr_glob][idx_quad], _normals[idx_cell][idx_nbr] );
                            break;
                        // second order solver
                        case 2:
                            // left status of interface
                            solL = _sol[idx_cell][idx_quad] +
                                   _limiter[idx_cell][idx_quad] *
                                       ( _solDx[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[idx_cell][0] ) +
                                         _solDy[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[idx_cell][1] ) );
                            solR = _sol[nbr_glob][idx_quad] +
                                   _limiter[nbr_glob][idx_quad] *
                                       ( _solDx[nbr_glob][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[nbr_glob][0] ) +
                                         _solDy[nbr_glob][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[nbr_glob][1] ) );
                            // flux evaluation
                            _solNew[idx_cell][idx_quad] += _g->Flux1D( _quadPoints[idx_quad], solL, solR, _normals[idx_cell][idx_nbr] );
                            break;
                            // higher order solver
                        default: ErrorMessages::Error( "Reconstruction order not supported.", CURRENT_FUNCTION ); break;
                    }
                }
            }
        }
    }
}

void SNSolver::FluxUpdatePseudo2D() {
    // Loop over all spatial cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        double solL;
        double solR;
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // Loop over all ordinates
        for( unsigned idx_quad = 0; idx_quad < _nq; ++idx_quad ) {
            // Reset temporary variable
            _solNew[idx_cell][idx_quad] = 0.0;
            // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_nbr = 0; idx_nbr < _neighbors[idx_cell].size(); ++idx_nbr ) {
                // store flux contribution on psiNew_sigmaS to save memory
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_nbr] == _nCells )
                    _solNew[idx_cell][idx_quad] +=
                        _g->Flux( _quadPoints[idx_quad], _sol[idx_cell][idx_quad], _sol[idx_cell][idx_quad], _normals[idx_cell][idx_nbr] );
                else {
                    unsigned int nbr_glob = _neighbors[idx_cell][idx_nbr];    // global idx of neighbor cell

                    switch( _reconsOrder ) {
                        // first order solver
                        case 1:
                            _solNew[idx_cell][idx_quad] +=
                                _g->Flux( _quadPoints[idx_quad], _sol[idx_cell][idx_quad], _sol[nbr_glob][idx_quad], _normals[idx_cell][idx_nbr] );
                            break;
                        // second order solver
                        case 2:
                            // left status of interface
                            solL = _sol[idx_cell][idx_quad] +
                                   _limiter[idx_cell][idx_quad] *
                                       ( _solDx[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[idx_cell][0] ) +
                                         _solDy[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[idx_cell][1] ) );
                            solR = _sol[nbr_glob][idx_quad] +
                                   _limiter[nbr_glob][idx_quad] *
                                       ( _solDx[nbr_glob][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[nbr_glob][0] ) +
                                         _solDy[nbr_glob][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[nbr_glob][1] ) );

                            // flux evaluation
                            _solNew[idx_cell][idx_quad] += _g->Flux( _quadPoints[idx_quad], solL, solR, _normals[idx_cell][idx_nbr] );
                            break;
                        default: ErrorMessages::Error( "Reconstruction order not supported.", CURRENT_FUNCTION ); break;
                    }
                }
            }
        }
    }
}

void SNSolver::FVMUpdate( unsigned idx_energy ) {
#pragma omp parallel for
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
