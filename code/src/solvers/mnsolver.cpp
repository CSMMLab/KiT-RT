#include "solvers/mnsolver.h"
#include "toolboxes/textprocessingtoolbox.h"
//#include <chrono>
#include <mpi.h>

MNSolver::MNSolver( Config* settings ) : Solver( settings ), _nMaxMomentsOrder( settings->GetMaxMomentDegree() ), _basis( _nMaxMomentsOrder ) {
    // Is this good (fast) code using a constructor list?

    _nTotalEntries = GlobalIndex( _nMaxMomentsOrder, int( _nMaxMomentsOrder ) ) + 1;
    _quadrature    = QuadratureBase::CreateQuadrature( _settings->GetQuadName(), settings->GetQuadOrder() );

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
    _moments = VectorVector( _nq, Vector( _nTotalEntries, 0.0 ) );
    ComputeMoments();
}

MNSolver::~MNSolver() {
    delete _quadrature;
    delete _entropy;
    delete _optimizer;
}

int MNSolver::GlobalIndex( int l, int k ) const {
    int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
    int prevIndicesThisLevel = k + l;    // number of previous indices in current level
    return numIndicesPrevLevel + prevIndicesThisLevel;
}

void MNSolver::ComputeMoments() {
    double x, y, z, w;

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        x                  = _quadrature->GetPoints()[idx_quad][0];
        y                  = _quadrature->GetPoints()[idx_quad][1];
        z                  = _quadrature->GetPoints()[idx_quad][2];
        _moments[idx_quad] = _basis.ComputeSphericalBasis( x, y, z );
    }

    Vector normalizationFactor( _moments[0].size(), 0.0 );

    // Normalize Moments

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        x = _quadrature->GetPoints()[idx_quad][0];
        y = _quadrature->GetPoints()[idx_quad][1];
        z = _quadrature->GetPoints()[idx_quad][2];
        w = _quadrature->GetWeights()[idx_quad];

        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            normalizationFactor[idx_sys] += w * _moments[idx_quad][idx_sys] * _moments[idx_quad][idx_sys];
        }
    }
    std::cout << "norm" << normalizationFactor << "\n";

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            _moments[idx_quad][idx_sys] /= sqrt( normalizationFactor[idx_sys] );
        }
    }
}

Vector MNSolver::ConstructFlux( unsigned idx_cell ) {

    // ---- Integration of Moment of flux ----
    double w, entropyL, entropyR, entropyFlux;
    // double x, y, z;
    // std::cout << "--------\n";
    // std::cout << " alphas curr cell" << _alpha[idx_cell] << " alphas neigh cell:\n";

    bool good = _alpha[idx_cell] == 0.0;

    int idx_good = 0;
    for( unsigned idx_neigh = 0; idx_neigh < _neighbors[idx_cell].size(); idx_neigh++ ) {
        // std::cout << _alpha[idx_neigh] << "\n";
        if( _alpha[idx_neigh] == 0.0 ) idx_good++;
    }
    if( idx_good == 2 && good ) {
        std::cout << "Candidate\n";
    }

    Vector flux( _nTotalEntries, 0.0 );
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {

        w = _quadrature->GetWeights()[idx_quad];

        entropyFlux = 0.0;    // Reset temorary flux

        entropyL = dot( _alpha[idx_cell], _moments[idx_quad] );
        // _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], _moments[idx_quad] ) );

        for( unsigned idx_neigh = 0; idx_neigh < _neighbors[idx_cell].size(); idx_neigh++ ) {
            // Store fluxes in psiNew, to save memory
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neigh] == _nCells )
                entropyR = entropyL;
            else {
                entropyR = dot( _alpha[idx_neigh], _moments[idx_quad] );
                //_entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_neigh], _moments[idx_quad] ) );
            }
            entropyFlux += _g->Flux( _quadrature->GetPoints()[idx_quad], entropyL, entropyR, _normals[idx_cell][idx_neigh] );
        }

        //   std::cout << "\n entropy Flux at cell [" << idx_cell << "] and quad [" << idx_quad << "]: " << entropyFlux << "\n";
        // integrate
        flux += _moments[idx_quad] * ( w * entropyFlux );
    }

    // std::cout << "\n integrate entropy Flux at cell [" << idx_cell << "] and : " << flux << "\n";

    return flux;
}

void MNSolver::Solve() {

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _psi;
    double dFlux        = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );

    double mass1 = 0;
    for( unsigned i = 0; i < _nCells; ++i ) {
        _solverOutput[i] = _psi[i][0];
        mass1 += _psi[i][0] * _areas[i];
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

            _alpha[idx_cell] = _psi[idx_cell];    // _optimizer->Solve( _psi[idx_cell] );

            // ------- Flux Computation Step ---------

            // Dirichlet Boundaries are finished now
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;

            psiNew[idx_cell] = ConstructFlux( idx_cell );

            // ------ FVM Step ------

            // NEED TO VECTORIZE
            for( unsigned idx_system = 0; idx_system < _nTotalEntries; idx_system++ ) {

                psiNew[idx_cell][idx_system] =
                    _psi[idx_cell][idx_system] - ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_system]; /* cell averaged flux */

                //- _dE * _psi[idx_cell][idx_system] *
                //      ( _sigmaA[idx_energy][idx_cell]                                    /* absorbtion influence */
                //        + _sigmaS[idx_energy][idx_cell] * _scatterMatDiag[idx_system] ); /* scattering influence */
            }
        }
        _psi = psiNew;

        // pseudo time iteration output
        double mass = 0.0;
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
            fluxNew[idx_cell]       = _psi[idx_cell][0];    // zeroth moment is raditation densitiy we are interested in
            _solverOutput[idx_cell] = _psi[idx_cell][0];
            mass += _psi[idx_cell][0] * _areas[idx_cell];
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
        flux[i] = _psi[i][0];
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
