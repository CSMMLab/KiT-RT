#include "solvers/mnsolver.h"
#include "toolboxes/textprocessingtoolbox.h"
#include <mpi.h>

MNSolver::MNSolver( Config* settings ) : Solver( settings ), _nMaxMomentsOrder( settings->GetMaxMomentDegree() ), _basis( _nMaxMomentsOrder ) {
    // Is this good (fast) code using a constructor list?

    _nTotalEntries = GlobalIndex( _nMaxMomentsOrder, int( _nMaxMomentsOrder ) ) + 1;
    _quadrature    = QuadratureBase::CreateQuadrature( settings->GetQuadName(), settings->GetQuadOrder() );
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

    // Initialize System Matrices
    _A = VectorVector( _nTotalEntries );

    _entropy = EntropyBase::Create( settings );
    // Initialize lagrange Multipliers
    _alpha = VectorVector( _nCells );
    for( unsigned idx_cells = 0; idx_cells < _nTotalEntries; idx_cells++ ) {
        _alpha[idx_cells] = Vector( _nTotalEntries, 0.0 );
    }
}

MNSolver::~MNSolver() { delete _quadrature; }

int MNSolver::GlobalIndex( int l, int k ) const {
    int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
    int prevIndicesThisLevel = k + l;    // number of previous indices in current level
    return numIndicesPrevLevel + prevIndicesThisLevel;
}

void MNSolver::ComputeSystemMatrices() {

    // Initialize entries of _A
    for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
        Vector entry( 2, 0.0 );
        _A[idx_sys] = entry;
    }

    // Compute Integrals _A[idx_system][i]=<Omega_i*m^(l,k)>
    Vector funcEval( _nTotalEntries, 0.0 );
    for( unsigned i = 0; i < _nq; i++ ) {
        double x = _quadrature->GetPoints()[i][0];
        double y = _quadrature->GetPoints()[i][1];
        double z = _quadrature->GetPoints()[i][2];
        double w = _quadrature->GetWeights()[i];
        funcEval = _basis.ComputeSphericalBasis( x, y, z );    // = m_l^k
        for( unsigned idx_system = 0; idx_system < _nTotalEntries; idx_system++ ) {
            _A[idx_system][0] += w * x * funcEval[idx_system];
            _A[idx_system][1] += w * y * funcEval[idx_system];
        }
    }
}

Vector MNSolver::ConstructFlux( unsigned idx_cell, unsigned idx_neigh ) {

    // ---- Integration of Moment of flux ----
    double x, y, z, w, entropy, velocity;
    Vector moment( _nTotalEntries, 0.0 );
    Vector flux( _nTotalEntries, 0.0 );
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        x      = _quadrature->GetPoints()[idx_quad][0];
        y      = _quadrature->GetPoints()[idx_quad][1];
        z      = _quadrature->GetPoints()[idx_quad][2];
        w      = _quadrature->GetWeights()[idx_quad];
        moment = _basis.ComputeSphericalBasis( x, y, z );

        // --- get upwinding direction (2D)
        // DOUBLE CHECK IF WE NEED DIMENSION SPLITTIN!!!

        velocity = _normals[idx_cell][idx_neigh][0] * x + _normals[idx_cell][idx_neigh][1] * y;

        // Perform upwinding
        ( velocity > 0 ) ? ( entropy = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], moment ) ) )
                         : ( entropy = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_neigh], moment ) ) );

        // integrate
        flux += moment * ( w * velocity * entropy );
    }

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
    unsigned idx_system = 0;
    double mass1        = 0;
    for( unsigned i = 0; i < _nCells; ++i ) {
        _solverOutput[i] = _psi[i][0];
        mass1 += _psi[i][0];
    }

    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "dFlux" );
    if( rank == 0 ) log->info( "{:03.8f}   {:01.5e} {:01.5e}", -1.0, dFlux, mass1 );

    // Loop over energies (pseudo-time of continuous slowing down approach)

    for( unsigned idx_energy = 0; idx_energy < _nEnergies; ++idx_energy ) {

        // Loop over the grid cells
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {

            // ------- Reconstruction Step -------

            _alpha[idx_cell] = _optimizer->Solve( _psi[idx_cell] );

            // ------- Compute Fluxes ---------

            // Loop over neighbors
            for( unsigned idx_neigh = 0; idx_neigh < _neighbors[idx_cell].size(); idx_neigh++ ) {
                // Store fluxes in psiNew, to save memory
                psiNew[idx_cell] += ConstructFlux( idx_cell, idx_neigh );
            }

            // ------ FVM Step ------
            // NEED TO VECTORIZE
            for( unsigned idx_system = 0; idx_system < _nTotalEntries; idx_system++ ) {

                psiNew[idx_cell][idx_system] = _psi[idx_cell][idx_system] -
                                               ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_system] /* cell averaged flux */
                                               - _dE * _psi[idx_cell][idx_system] *
                                                     ( _sigmaA[idx_energy][idx_cell]                                    /* absorbtion influence */
                                                       + _sigmaS[idx_energy][idx_cell] * _scatterMatDiag[idx_system] ); /* scattering influence */
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

        Save( idx_energy );
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e}", _energies[idx_energy], dFlux );
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
