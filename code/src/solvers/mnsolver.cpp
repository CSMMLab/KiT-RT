#include "solvers/mnsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "entropies/entropybase.h"
#include "fluxes/numericalflux.h"
#include "optimizers/optimizerbase.h"
#include "quadratures/quadraturebase.h"
#include "solvers/sphericalharmonics.h"
#include "toolboxes/textprocessingtoolbox.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

#include <fstream>
//#include <chrono>

MNSolver::MNSolver( Config* settings ) : Solver( settings ) {

    // Is this good (fast) code using a constructor list?
    _nMaxMomentsOrder = settings->GetMaxMomentDegree();
    _nTotalEntries    = GlobalIndex( _nMaxMomentsOrder, int( _nMaxMomentsOrder ) ) + 1;

    // build quadrature object and store quadrature points and weights
    _quadPoints       = _quadrature->GetPoints();
    _weights          = _quadrature->GetWeights();
    _nq               = _quadrature->GetNq();
    _quadPointsSphere = _quadrature->GetPointsSphere();
    _settings->SetNQuadPoints( _nq );

    // transform sigmaT and sigmaS in sigmaA.
    _sigmaA = VectorVector( _nEnergies, Vector( _nCells, 0 ) );    // Get rid of this extra vektor!

    for( unsigned n = 0; n < _nEnergies; n++ ) {
        for( unsigned j = 0; j < _nCells; j++ ) {
            _sigmaA[n][j] = 0;    //_sigmaT[n][j] - _sigmaS[n][j];
            _sigmaS[n][j] = 1;
        }
    }

    // Initialize Scatter Matrix
    _scatterMatDiag    = Vector( _nTotalEntries, 1.0 );
    _scatterMatDiag[0] = 0.0;    // First entry is zero by construction.

    // Initialize Entropy
    _entropy = EntropyBase::Create( _settings );

    // Initialize Optimizer
    _optimizer = OptimizerBase::Create( _settings );

    // Initialize lagrange Multiplier
    _alpha = VectorVector( _nCells, Vector( _nTotalEntries, 0.0 ) );

    // Initialize and Pre-Compute Moments at quadrature points
    _basis = new SphericalHarmonics( _nMaxMomentsOrder );

    _moments = VectorVector( _nq, Vector( _nTotalEntries, 0.0 ) );
    ComputeMoments();

    // Solver output
    _outputFields = std::vector( _nTotalEntries, std::vector( _nCells, 0.0 ) );
}

MNSolver::~MNSolver() {
    delete _entropy;
    delete _optimizer;
    delete _basis;
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

void MNSolver::Solve() {

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _sol;

    double mass = 0;

    Save( -1 );

    if( rank == 0 ) log->info( "{:10}  {:10}", "t", "mass" );

    // if( rank == 0 ) log->info( "{:03.8f}   {:01.5e} {:01.5e}", -1.0, dFlux, mass1 ); Should be deleted

    // Loop over energies (pseudo-time of continuous slowing down approach)
    for( unsigned idx_energy = 0; idx_energy < _nEnergies; idx_energy++ ) {

        // ------- Reconstruction Step ----------------

        _optimizer->SolveMultiCell( _alpha, _sol, _moments );

        // Loop over the grid cells
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {

            // ------- Relizablity Reconstruction Step ----

            ComputeRealizableSolution( idx_cell );

            // ------- Flux Computation Step --------------

            // Dirichlet Boundaries are finished now
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;

            psiNew[idx_cell] = ConstructFlux( idx_cell );

            // ------ Finite Volume Update Step ------

            // NEED TO VECTORIZE
            for( unsigned idx_system = 0; idx_system < _nTotalEntries; idx_system++ ) {

                psiNew[idx_cell][idx_system] = _sol[idx_cell][idx_system] -
                                               ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_system] /* cell averaged flux */
                                               - _dE * _sol[idx_cell][idx_system] *
                                                     ( _sigmaA[idx_energy][idx_cell]                                    /* absorbtion influence */
                                                       + _sigmaS[idx_energy][idx_cell] * _scatterMatDiag[idx_system] ); /* scattering influence */
            }
        }

        // Update Solution
        _sol = psiNew;

        // --- VTK and CSV Output ---
        mass = WriteOutputFields();
        Save( idx_energy );
        WriteNNTrainingData( idx_energy );

        // --- Screen Output ---
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e}", _energies[idx_energy], mass );
    }
}

double MNSolver::WriteOutputFields() {
    double mass = 0.0;

    for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
            mass += _sol[idx_cell][0] * _areas[idx_cell];
            _outputFields[idx_sys][idx_cell] = _sol[idx_cell][idx_sys];
        }
    }

    return mass;
}

void MNSolver::Save() const { Save( _nEnergies - 1 ); }

void MNSolver::Save( int currEnergy ) const {
    std::vector<std::string> fieldNames( _nTotalEntries, "flux" );
    for( int idx_l = 0; idx_l <= (int)_nMaxMomentsOrder; idx_l++ ) {
        for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
            fieldNames[GlobalIndex( idx_l, idx_k )] = std::string( "u_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
        }
    }
    std::vector<std::vector<std::vector<double>>> results{ _outputFields };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), results, fieldNames, _mesh );
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
