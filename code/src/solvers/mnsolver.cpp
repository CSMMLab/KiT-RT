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
    }
    return flux;
}

void MNSolver::Solve() {

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _sol;
    double dFlux        = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );

    double mass1 = 0;
    for( unsigned i = 0; i < _nCells; ++i ) {
        _solverOutput[i] = _sol[i][0];
        mass1 += _sol[i][0] * _areas[i];
    }

    dFlux   = blaze::l2Norm( fluxNew - fluxOld );
    fluxOld = fluxNew;

    Save( -1 );

    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "dFlux" );
    if( rank == 0 ) log->info( "{:03.8f}   {:01.5e} {:01.5e}", -1.0, dFlux, mass1 );

    // Time measurement
    // auto start = chrono::steady_clock::now();
    // auto end   = chrono::steady_clock::now();

    // Loop over energies (pseudo-time of continuous slowing down approach)

    for( unsigned idx_energy = 0; idx_energy < _nEnergies; idx_energy++ ) {

        // Loop over the grid cells
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {

            // ------- Reconstruction Step -------

            // _alpha[idx_cell] = _sol[idx_cell];
            _optimizer->Solve( _alpha[idx_cell], _sol[idx_cell], _moments );

            // ------- Flux Computation Step ---------

            // Dirichlet Boundaries are finished now
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;

            psiNew[idx_cell] = ConstructFlux( idx_cell );

            // ------ FVM Step ------

            // NEED TO VECTORIZE
            for( unsigned idx_system = 0; idx_system < _nTotalEntries; idx_system++ ) {

                psiNew[idx_cell][idx_system] = _sol[idx_cell][idx_system] -
                                               ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_system] /* cell averaged flux */
                                               - _dE * _sol[idx_cell][idx_system] *
                                                     ( _sigmaA[idx_energy][idx_cell]                                    /* absorbtion influence */
                                                       + _sigmaS[idx_energy][idx_cell] * _scatterMatDiag[idx_system] ); /* scattering influence */
            }
        }
        _sol = psiNew;

        // pseudo time iteration output
        double mass = 0.0;
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
            fluxNew[idx_cell]       = _sol[idx_cell][0];    // zeroth moment is raditation densitiy we are interested in
            _solverOutput[idx_cell] = _sol[idx_cell][0];
            mass += _sol[idx_cell][0] * _areas[idx_cell];
        }

        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e} {:01.5e}", _energies[idx_energy], dFlux, mass );
        Save( idx_energy );
    }
}

void MNSolver::Save() const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<double> flux;
    flux.resize( _nCells );

    for( unsigned i = 0; i < _nCells; ++i ) {
        flux[i] = _sol[i][0];
    }
    std::vector<std::vector<double>> scalarField( 1, flux );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile(), results, fieldNames, _mesh );
}

void MNSolver::Save( int currEnergy ) const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<std::vector<double>> scalarField( 1, _solverOutput );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), results, fieldNames, _mesh );
}