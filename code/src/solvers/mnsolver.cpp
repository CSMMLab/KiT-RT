#include "solvers/mnsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "entropies/entropybase.h"
#include "fluxes/numericalflux.h"
#include "optimizers/optimizerbase.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"
#include "solvers/sphericalharmonics.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/textprocessingtoolbox.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

#include <fstream>
//#include <chrono>

MNSolver::MNSolver( Config* settings ) : Solver( settings ) {

    // Is this good (fast) code using a constructor list?
    _LMaxDegree    = settings->GetMaxMomentDegree();
    _nTotalEntries = GlobalIndex( _LMaxDegree, int( _LMaxDegree ) ) + 1;

    // build quadrature object and store quadrature points and weights
    _quadPoints       = _quadrature->GetPoints();
    _weights          = _quadrature->GetWeights();
    _nq               = _quadrature->GetNq();
    _quadPointsSphere = _quadrature->GetPointsSphere();
    _settings->SetNQuadPoints( _nq );

    // Initialize Scatter Matrix --
    _scatterMatDiag = Vector( _nTotalEntries, 0.0 );
    ComputeScatterMatrix();

    // Initialize Entropy
    _entropy = EntropyBase::Create( _settings );

    // Initialize Optimizer
    _optimizer = OptimizerBase::Create( _settings );

    // Initialize lagrange Multiplier
    _alpha = VectorVector( _nCells, Vector( _nTotalEntries, 0.0 ) );

    // Initialize and Pre-Compute Moments at quadrature points
    _basis = new SphericalHarmonics( _LMaxDegree );

    _moments = VectorVector( _nq, Vector( _nTotalEntries, 0.0 ) );
    ComputeMoments();

    // Solver output
    PrepareVolumeOutput();
}

MNSolver::~MNSolver() {
    delete _entropy;
    delete _optimizer;
    delete _basis;
}

void MNSolver::ComputeScatterMatrix() {

    // --- Isotropic ---
    _scatterMatDiag[0] = -1.0;
    for( unsigned idx_diag = 1; idx_diag < _nTotalEntries; idx_diag++ ) {
        _scatterMatDiag[idx_diag] = 0.0;
    }
}

int MNSolver::GlobalIndex( int l, int k ) const {
    int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
    int prevIndicesThisLevel = k + l;    // number of previous indices in current level
    return numIndicesPrevLevel + prevIndicesThisLevel;
}

void MNSolver::ComputeMoments() {
    double my, phi;

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my  = _quadPointsSphere[idx_quad][0];
        phi = _quadPointsSphere[idx_quad][1];

        _moments[idx_quad] = _basis->ComputeSphericalBasis( my, phi );
    }
}

Vector MNSolver::ConstructFlux( unsigned idx_cell ) {

    // ---- Integration of Moment of flux ----
    double entropyL, entropyR, entropyFlux;

    Vector flux( _nTotalEntries, 0.0 );

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {

        entropyFlux = 0.0;    // Reset temorary flux

        entropyL = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], _moments[idx_quad] ) );

        for( unsigned idx_neigh = 0; idx_neigh < _neighbors[idx_cell].size(); idx_neigh++ ) {
            // Store fluxes in psiNew, to save memory
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neigh] == _nCells )
                entropyR = entropyL;
            else {
                entropyR = _entropy->EntropyPrimeDual( blaze::dot( _alpha[_neighbors[idx_cell][idx_neigh]], _moments[idx_quad] ) );
            }
            entropyFlux += _g->Flux( _quadPoints[idx_quad], entropyL, entropyR, _normals[idx_cell][idx_neigh] );
        }
        flux += _moments[idx_quad] * ( _weights[idx_quad] * entropyFlux );

        // ------- Relizablity Reconstruction Step ----
    }
    return flux;
}

void MNSolver::ComputeRealizableSolution( unsigned idx_cell ) {
    double entropyReconstruction = 0.0;
    _sol[idx_cell]               = 0;
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        // Make entropyReconstruction a member vector, s.t. it does not have to be re-evaluated in ConstructFlux
        entropyReconstruction = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], _moments[idx_quad] ) );
        _sol[idx_cell] += _moments[idx_quad] * ( _weights[idx_quad] * entropyReconstruction );
    }
}

void MNSolver::IterPreprocessing() {

    // ------- Reconstruction Step ----------------

    _optimizer->SolveMultiCell( _alpha, _sol, _moments );

    // ------- Relizablity Reconstruction Step ----
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        ComputeRealizableSolution( idx_cell );
    }
}

void MNSolver::IterPostprocessing() {
    // --- Update Solution ---
    _sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();
}

void MNSolver::ComputeRadFlux() {
    double firstMomentScaleFactor = sqrt( 4 * M_PI );
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        _fluxNew[idx_cell] = _sol[idx_cell][0] * firstMomentScaleFactor;
    }
}

void MNSolver::FluxUpdate() {
    // Loop over the grid cells
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        // Dirichlet Boundaries stay
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        _solNew[idx_cell] = ConstructFlux( idx_cell );
    }
}

void MNSolver::FVMUpdate( unsigned idx_energy ) {
    // Loop over the grid cells
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        // Dirichlet Boundaries stay
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;

        for( unsigned idx_system = 0; idx_system < _nTotalEntries; idx_system++ ) {

            _solNew[idx_cell][idx_system] = _sol[idx_cell][idx_system] -
                                            ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_system] /* cell averaged flux */
                                            - _dE * _sol[idx_cell][idx_system] *
                                                  ( _sigmaT[idx_energy][idx_cell]                                    /* absorbtion influence */
                                                    + _sigmaS[idx_energy][idx_cell] * _scatterMatDiag[idx_system] ); /* scattering influence */
        }

        _solNew[idx_cell][0] += _dE * _Q[0][idx_cell][0];
    }
}

void MNSolver::PrepareVolumeOutput() {
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

            case MOMENTS:
                // As many entries as there are moments in the system
                _outputFields[idx_group].resize( _nTotalEntries );
                _outputFieldNames[idx_group].resize( _nTotalEntries );

                for( int idx_l = 0; idx_l <= (int)_LMaxDegree; idx_l++ ) {
                    for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                        _outputFields[idx_group][GlobalIndex( idx_l, idx_k )].resize( _nCells );

                        _outputFieldNames[idx_group][GlobalIndex( idx_l, idx_k )] =
                            std::string( "u_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                    }
                }
                break;

            case DUAL_MOMENTS:
                // As many entries as there are moments in the system
                _outputFields[idx_group].resize( _nTotalEntries );
                _outputFieldNames[idx_group].resize( _nTotalEntries );

                for( int idx_l = 0; idx_l <= (int)_LMaxDegree; idx_l++ ) {
                    for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                        _outputFields[idx_group][GlobalIndex( idx_l, idx_k )].resize( _nCells );

                        _outputFieldNames[idx_group][GlobalIndex( idx_l, idx_k )] =
                            std::string( "alpha_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                    }
                }
                break;

            case ANALYTIC:
                // one entry per cell
                _outputFields[idx_group].resize( 1 );
                _outputFieldNames[idx_group].resize( 1 );
                _outputFields[idx_group][0].resize( _nCells );
                _outputFieldNames[idx_group][0] = std::string( "analytic radiation flux density" );
                break;

            default: ErrorMessages::Error( "Volume Output Group not defined for MN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

void MNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();

    // Check if volume output fields are written to file this iteration
    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_pseudoTime % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_pseudoTime == _nEnergies - 1 ) /* need sol at last iteration */ ) {

        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _fluxNew[idx_cell];
                    }
                    break;
                case MOMENTS:
                    for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
                        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                            _outputFields[idx_group][idx_sys][idx_cell] = _sol[idx_cell][idx_sys];
                        }
                    }
                    break;
                case DUAL_MOMENTS:
                    for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
                        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                            _outputFields[idx_group][idx_sys][idx_cell] = _alpha[idx_cell][idx_sys];
                        }
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

                default: ErrorMessages::Error( "Volume Output Group not defined for MN Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}

void MNSolver::WriteNNTrainingData( unsigned idx_pseudoTime ) {
    std::string filename = "trainNN.csv";
    std::ofstream myfile;
    myfile.open( filename, std::ofstream::app );

    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        myfile << 0 << ", " << _nTotalEntries << "," << idx_pseudoTime;
        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            myfile << "," << _sol[idx_cell][idx_sys];
        }
        myfile << " \n" << 1 << ", " << _nTotalEntries << "," << idx_pseudoTime;
        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            myfile << "," << _alpha[idx_cell][idx_sys];
        }
        myfile << "\n";
    }
    myfile.close();
}
