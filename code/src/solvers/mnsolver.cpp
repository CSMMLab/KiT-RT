#include "solvers/mnsolver.h"

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

    // Fill System Matrices
    ComputeSystemMatrices();
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
    std::vector<double> funcEval( _nTotalEntries, 0.0 );
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

        // ------- Reconstruction Step -------
        // Todo
        // ------- Advection Step ------------

        // Loop over all spatial cells
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;    // Dirichlet cells stay at IC, farfield assumption

            // Loop over all equations of the system.
            for( int idx_lDegree = 0; idx_lDegree <= int( _nMaxMomentsOrder ); idx_lDegree++ ) {
                for( int idx_kOrder = -idx_lDegree; idx_kOrder <= idx_lDegree; idx_kOrder++ ) {
                    idx_system                   = unsigned( GlobalIndex( idx_lDegree, idx_kOrder ) );
                    psiNew[idx_cell][idx_system] = 0.0;

                    // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                    for( unsigned l = 0; l < _neighbors[idx_cell].size(); ++l ) {
                        // store flux contribution on psiNew_sigmaS to save memory
                        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][l] == _nCells )
                            psiNew[idx_cell][idx_system] +=
                                _g->Flux( _A[idx_system], _psi[idx_cell][idx_system], _psi[idx_cell][idx_system], _normals[idx_cell][l] );
                        else
                            psiNew[idx_cell][idx_system] += _g->Flux(
                                _A[idx_system], _psi[idx_cell][idx_system], _psi[_neighbors[idx_cell][l]][idx_system], _normals[idx_cell][l] );
                    }
                    // time update angular flux with numerical flux and total scattering cross section
                    psiNew[idx_cell][idx_system] = _psi[idx_cell][idx_system] -
                                                   ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_system] /* cell averaged flux */
                                                   - _dE * _psi[idx_cell][idx_system] *
                                                         ( _sigmaA[idx_energy][idx_cell]                                    /* absorbtion influence */
                                                           + _sigmaS[idx_energy][idx_cell] * _scatterMatDiag[idx_system] ); /* scattering influence */
                }
            }

            // TODO: figure out a more elegant way
            // add external source contribution
            // if( _Q.size() == 1u ) {                   // constant source for all energies
            //    if( _Q[0][idx_cell].size() == 1u )    // isotropic source
            //        psiNew[idx_cell] += _dE * _Q[0][idx_cell][0];
            //    else
            //        psiNew[idx_cell] += _dE * _Q[0][idx_cell];
            //}
            // else {
            //    if( _Q[0][idx_cell].size() == 1u )    // isotropic source
            //        psiNew[idx_cell] += _dE * _Q[idx_energy][idx_cell][0];
            //    else
            //        psiNew[idx_cell] += _dE * _Q[idx_energy][idx_cell];
            //}
        }
        _psi        = psiNew;
        double mass = 0.0;
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
            fluxNew[idx_cell]       = _psi[idx_cell][0];    // zeroth moment is raditation densitiy we are interested in
            _solverOutput[idx_cell] = _psi[idx_cell][0];
            mass += _psi[idx_cell][0] * _areas[idx_cell];
        }
        // Save( n );
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e}", _energies[idx_energy], dFlux );
    }
}
