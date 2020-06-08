
#include "pnsolver.h"
#include "../ext/blaze/blaze/math/blas/cblas/gemm.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/textprocessingtoolbox.h"
#include <mpi.h>

PNSolver::PNSolver( Config* settings ) : Solver( settings ) {
    _nTotalEntries = GlobalIndex( _nq, int( _nq ) ) + 1;

    // transform sigmaT and sigmaS in sigmaA.
    _sigmaA = VectorVector( _nEnergies, Vector( _nCells, 0 ) );    // Get rid of this extra vektor!
    for( unsigned n = 0; n < _nEnergies; n++ ) {
        for( unsigned j = 0; j < _nCells; j++ ) {
            _sigmaA[n][j] = 0;    //_sigmaT[n][j] - _sigmaS[n][j];
            _sigmaS[n][j] = 1;
        }
    }

    // Initialize System Matrices
    _Ax = SymMatrix( _nTotalEntries );
    _Ay = SymMatrix( _nTotalEntries );
    _Az = SymMatrix( _nTotalEntries );

    _AxPlus  = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AxMinus = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AxAbs   = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AyPlus  = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AyMinus = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AyAbs   = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AzPlus  = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AzMinus = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AzAbs   = Matrix( _nTotalEntries, _nTotalEntries, 0 );

    // Initialize Scatter Matrix
    _scatterMatDiag = Vector( _nTotalEntries, 0 );

    // Fill System Matrices
    ComputeSystemMatrices();

    std::cout << "System Matrix Set UP!" << std::endl;
    // Compute Decomposition in positive and negative (eigenvalue) parts of flux jacobians
    ComputeFluxComponents();
    std::cout << "_Ax :\n" << _Ax << "\n ";    // _AxP \n" << _AxPlus << "\n _AxM \n" << _AxMinus << "\n";
    std::cout << "_Ay :\n" << _Ay << "\n ";    //_AyP \n" << _AyPlus << "\n _AyM \n" << _AyMinus << "\n";
    std::cout << "_Az :\n" << _Az << "\n ";    //_AzP \n" << _AzPlus << "\n _AzM \n" << _AzMinus << "\n";

    std::cout << "_AxA :\n" << _AxAbs << "\n ";    // _AxP \n" << _AxPlus << "\n _AxM \n" << _AxMinus << "\n";
    std::cout << "_AyA :\n" << _AyAbs << "\n ";    //_AyP \n" << _AyPlus << "\n _AyM \n" << _AyMinus << "\n";
    std::cout << "_AzA :\n" << _AzAbs << "\n ";

    std::cout << "_AxA :\n" << _AxPlus - _AxMinus << "\n ";    // _AxP \n" << _AxPlus << "\n _AxM \n" << _AxMinus << "\n";
    std::cout << "_AyA :\n" << _AyPlus - _AyMinus << "\n ";    //_AyP \n" << _AyPlus << "\n _AyM \n" << _AyMinus << "\n";
    std::cout << "_AzA :\n" << _AzPlus - _AzMinus << "\n ";

    // Compute diagonal of the scatter matrix (it's a diagonal matrix)
    ComputeScatterMatrix();

    std::cout << "scatterMatrix : " << _scatterMatDiag << "\n";
}

double PNSolver::CTilde( int l, int k ) const {
    if( k < 0 ) return 0.0;
    if( k == 0 )
        return std::sqrt( 2 ) * CParam( l, k );
    else
        return CParam( l, k );
}

double PNSolver::DTilde( int l, int k ) const {
    if( k < 0 ) return 0.0;
    if( k == 0 )
        return std::sqrt( 2 ) * DParam( l, k );
    else
        return DParam( l, k );
}

double PNSolver::ETilde( int l, int k ) const {
    if( k == 1 )
        return std::sqrt( 2 ) * EParam( l, k );
    else
        return EParam( l, k );
}

double PNSolver::FTilde( int l, int k ) const {
    if( k == 1 )
        return std::sqrt( 2 ) * FParam( l, k );
    else
        return FParam( l, k );
}

double PNSolver::AParam( int l, int k ) const {
    return std::sqrt( double( ( l - k + 1 ) * ( l + k + 1 ) ) / double( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) );
}

double PNSolver::BParam( int l, int k ) const { return std::sqrt( double( ( l - k ) * ( l + k ) ) / double( ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ) ); }

double PNSolver::CParam( int l, int k ) const {
    return std::sqrt( double( ( l + k + 1 ) * ( l + k + 2 ) ) / double( ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ) );
}

double PNSolver::DParam( int l, int k ) const {
    return std::sqrt( double( ( l - k ) * ( l - k - 1 ) ) / double( ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ) );
}

double PNSolver::EParam( int l, int k ) const {
    return std::sqrt( double( ( l - k + 1 ) * ( l - k + 2 ) ) / double( ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ) );
}

double PNSolver::FParam( int l, int k ) const { return std::sqrt( double( ( l + k ) * ( l + k - 1 ) ) / double( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ); }

int PNSolver::kPlus( int k ) const { return k + Sgn( k ); }

int PNSolver::kMinus( int k ) const { return k - Sgn( k ); }

int PNSolver::GlobalIndex( int l, int k ) const {
    int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
    int prevIndicesThisLevel = k + l;    // number of previous indices in current level
    return numIndicesPrevLevel + prevIndicesThisLevel;
}

bool PNSolver::CheckIndex( int l, int k ) const {
    if( l >= 0 && l <= int( _nq ) ) {
        if( k >= -l && k <= l ) return true;
    }
    return false;
}

int PNSolver::Sgn( int k ) const {
    if( k >= 0 )
        return 1;
    else
        return -1;
}

void PNSolver::ComputeSystemMatrices() {
    int idx_col      = 0;
    unsigned idx_row = 0;

    // loop over columns of A
    for( int idx_lOrder = 0; idx_lOrder <= int( _nq ); idx_lOrder++ ) {                  // index of legendre polynom
        for( int idx_kOrder = -idx_lOrder; idx_kOrder <= idx_lOrder; idx_kOrder++ ) {    // second index of legendre function
            idx_row = unsigned( GlobalIndex( idx_lOrder, idx_kOrder ) );

            // flux matrix in direction x
            {
                if( idx_kOrder != -1 ) {

                    if( CheckIndex( idx_lOrder - 1, kMinus( idx_kOrder ) ) ) {
                        idx_col                             = GlobalIndex( idx_lOrder - 1, kMinus( idx_kOrder ) );
                        _Ax( idx_row, unsigned( idx_col ) ) = 0.5 * CTilde( idx_lOrder - 1, std::abs( idx_kOrder ) - 1 );
                    }

                    if( CheckIndex( idx_lOrder + 1, kMinus( idx_kOrder ) ) ) {
                        idx_col                             = GlobalIndex( idx_lOrder + 1, kMinus( idx_kOrder ) );
                        _Ax( idx_row, unsigned( idx_col ) ) = -0.5 * DTilde( idx_lOrder + 1, std::abs( idx_kOrder ) - 1 );
                    }
                }

                if( CheckIndex( idx_lOrder - 1, kPlus( idx_kOrder ) ) ) {
                    idx_col                             = GlobalIndex( idx_lOrder - 1, kPlus( idx_kOrder ) );
                    _Ax( idx_row, unsigned( idx_col ) ) = -0.5 * ETilde( idx_lOrder - 1, std::abs( idx_kOrder ) + 1 );
                }

                if( CheckIndex( idx_lOrder + 1, kPlus( idx_kOrder ) ) ) {
                    idx_col                             = GlobalIndex( idx_lOrder + 1, kPlus( idx_kOrder ) );
                    _Ax( idx_row, unsigned( idx_col ) ) = 0.5 * FTilde( idx_lOrder + 1, std::abs( idx_kOrder ) + 1 );
                }
            }

            // flux matrix in direction y
            {
                if( idx_kOrder != 1 ) {
                    if( CheckIndex( idx_lOrder - 1, -kMinus( idx_kOrder ) ) ) {
                        idx_col                             = GlobalIndex( idx_lOrder - 1, -kMinus( idx_kOrder ) );
                        _Ay( idx_row, unsigned( idx_col ) ) = -0.5 * Sgn( idx_kOrder ) * CTilde( idx_lOrder - 1, std::abs( idx_kOrder ) - 1 );
                    }

                    if( CheckIndex( idx_lOrder + 1, -kMinus( idx_kOrder ) ) ) {
                        idx_col                             = GlobalIndex( idx_lOrder + 1, -kMinus( idx_kOrder ) );
                        _Ay( idx_row, unsigned( idx_col ) ) = 0.5 * Sgn( idx_kOrder ) * DTilde( idx_lOrder + 1, std::abs( idx_kOrder ) - 1 );
                    }
                }

                if( CheckIndex( idx_lOrder - 1, -kPlus( idx_kOrder ) ) ) {
                    idx_col                             = GlobalIndex( idx_lOrder - 1, -kPlus( idx_kOrder ) );
                    _Ay( idx_row, unsigned( idx_col ) ) = -0.5 * Sgn( idx_kOrder ) * ETilde( idx_lOrder - 1, std::abs( idx_kOrder ) + 1 );
                }

                if( CheckIndex( idx_lOrder + 1, -kPlus( idx_kOrder ) ) ) {
                    idx_col                             = GlobalIndex( idx_lOrder + 1, -kPlus( idx_kOrder ) );
                    _Ay( idx_row, unsigned( idx_col ) ) = 0.5 * Sgn( idx_kOrder ) * FTilde( idx_lOrder + 1, std::abs( idx_kOrder ) + 1 );
                }
            }

            // flux matrix in direction z
            {
                if( CheckIndex( idx_lOrder - 1, idx_kOrder ) ) {
                    idx_col                             = GlobalIndex( idx_lOrder - 1, idx_kOrder );
                    _Az( idx_row, unsigned( idx_col ) ) = AParam( idx_lOrder - 1, idx_kOrder );
                }

                if( CheckIndex( idx_lOrder + 1, idx_kOrder ) ) {
                    idx_col                             = GlobalIndex( idx_lOrder + 1, idx_kOrder );
                    _Az( idx_row, unsigned( idx_col ) ) = BParam( idx_lOrder + 1, idx_kOrder );
                }
            }
        }
    }
}

void PNSolver::ComputeFluxComponents() {
    Vector eigenValues( _nTotalEntries, 0 );
    Vector eigenValuesX( _nTotalEntries, 0 );
    Vector eigenValuesY( _nTotalEntries, 0 );

    MatrixCol eigenVectors( _nTotalEntries, _nTotalEntries, 0 );    // ColumnMatrix for _AxPlus * eigenVectors Multiplication via SIMD
    // --- For x Direction ---
    {
        blaze::eigen( _Ax, eigenValues, eigenVectors );    // Compute Eigenvalues and Eigenvectors

        // Compute Flux Matrices A+ and A-
        for( unsigned idx_ij = 0; idx_ij < _nTotalEntries; idx_ij++ ) {
            if( eigenValues[idx_ij] >= 0 ) {
                _AxPlus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // positive part of Diagonal Matrix stored in _AxPlus
                _AxAbs( idx_ij, idx_ij )  = eigenValues[idx_ij];
            }
            else {
                _AxMinus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // negative part of Diagonal Matrix stored in _AxMinus
                _AxAbs( idx_ij, idx_ij )   = -eigenValues[idx_ij];
            }
        }

        _AxPlus  = eigenVectors * _AxPlus;    // col*row minimum performance
        _AxMinus = eigenVectors * _AxMinus;
        _AxAbs   = eigenVectors * _AxAbs;
        blaze::transpose( eigenVectors );
        _AxPlus  = _AxPlus * eigenVectors;    // row*col maximum performance
        _AxMinus = _AxMinus * eigenVectors;
        _AxAbs   = _AxAbs * eigenVectors;

        eigenValuesX = eigenValues;
    }
    // --- For y Direction -------
    {
        blaze::eigen( _Ay, eigenValues, eigenVectors );    // Compute Eigenvalues and Eigenvectors

        // Compute Flux Matrices A+ and A-
        for( unsigned idx_ij = 0; idx_ij < _nTotalEntries; idx_ij++ ) {
            if( eigenValues[idx_ij] >= 0 ) {
                _AyPlus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // positive part of Diagonal Matrix stored in _AxPlus
                _AyAbs( idx_ij, idx_ij )  = eigenValues[idx_ij];
            }
            else {
                _AyMinus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // negative part of Diagonal Matrix stored in _AxMinus
                _AyAbs( idx_ij, idx_ij )   = -eigenValues[idx_ij];
            }
        }

        _AyPlus  = eigenVectors * _AyPlus;
        _AyMinus = eigenVectors * _AyMinus;
        _AyAbs   = eigenVectors * _AyAbs;
        blaze::transpose( eigenVectors );
        _AyPlus  = _AyPlus * eigenVectors;
        _AyMinus = _AyMinus * eigenVectors;
        _AyAbs   = _AyAbs * eigenVectors;

        eigenValuesY = eigenValues;
    }
    // --- For z Direction -------
    {
        blaze::eigen( _Az, eigenValues, eigenVectors );    // Compute Eigenvalues and Eigenvectors

        // Compute Flux Matrices A+ and A-
        for( unsigned idx_ij = 0; idx_ij < _nTotalEntries; idx_ij++ ) {
            if( eigenValues[idx_ij] >= 0 ) {
                _AzPlus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // positive part of Diagonal Matrix stored in _AxPlus
                _AzAbs( idx_ij, idx_ij )  = eigenValues[idx_ij];
            }
            else {
                _AzMinus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // negative part of Diagonal Matrix stored in _AxMinus
                _AzAbs( idx_ij, idx_ij )   = -eigenValues[idx_ij];
            }
        }

        _AzPlus  = eigenVectors * _AzPlus;
        _AzMinus = eigenVectors * _AzMinus;
        _AzAbs   = eigenVectors * _AzAbs;
        blaze::transpose( eigenVectors );
        _AzPlus  = _AzPlus * eigenVectors;
        _AzMinus = _AzMinus * eigenVectors;
        _AzAbs   = _AzAbs * eigenVectors;
    }

    // Compute Spectral Radius
    std::cout << "Eigenvalues x direction " << eigenValuesX << "\n";
    std::cout << "Eigenvalues y direction " << eigenValuesY << "\n";
    std::cout << "Eigenvalues z direction " << eigenValues << "\n";

    std::cout << "Spectral Radius X " << blaze::max( blaze::abs( eigenValuesX ) ) << "\n";
    std::cout << "Spectral Radius Y " << blaze::max( blaze::abs( eigenValuesY ) ) << "\n";
    std::cout << "Spectral Radius Z " << blaze::max( blaze::abs( eigenValues ) ) << "\n";

    std::cout << "Spectral Radius combined " << blaze::max( blaze::abs( eigenValues + eigenValuesX + eigenValuesY ) ) << "\n";
}

void PNSolver::ComputeScatterMatrix() {

    // right now independent of spatial coordinates
    // G[i][i] = 1- g_l
    // g_l = 2*pi*\int_-1 ^1 k(x,my)P_l(my)dmy
    // use midpoint rule.
    // Here is room for optimization!

    // Only isotropic scattering: k = 1/(2pi) in 2d !!!

    unsigned nSteps     = 500000;
    double integralSum  = 0;
    unsigned idx_global = 0;

    for( int idx_lOrder = 0; idx_lOrder <= int( _nq ); idx_lOrder++ ) {
        for( int idx_kOrder = -idx_lOrder; idx_kOrder <= idx_lOrder; idx_kOrder++ ) {
            idx_global                  = unsigned( GlobalIndex( idx_lOrder, idx_kOrder ) );
            _scatterMatDiag[idx_global] = 0;

            // compute value of diagonal
            integralSum = 0.5 * ( Legendre( -1.0, idx_lOrder ) + Legendre( 1.0, idx_lOrder ) );    // Boundary terms

            for( unsigned i = 1; i < nSteps - 1; i++ ) {
                integralSum += Legendre( -1.0 + double( i ) * 2 / double( nSteps ), idx_lOrder );
            }
            integralSum *= 2 / double( nSteps );

            // isotropic scattering and prefactor of eigenvalue
            integralSum *= 0.5;    // 2pi/4pi
            _scatterMatDiag[idx_global] = 1 - integralSum;
        }
    }
}

double PNSolver::Legendre( double x, int l ) {
    // Pre computed low order polynomials for faster computation
    switch( l ) {
        case 0: return 1;
        case 1: return x;
        case 2:    // 0.5*(3*x*x - 1)
            return 1.5 * x * x - 0.5;
        case 3:    // 0.5* (5*x*x*x -3 *x)
            return 2.5 * x * x * x - 1.5 * x;
        case 4:    // 1/8*(35x^4-30x^2 + 3)
            return 4.375 * x * x * x * x - 3.75 * x * x + 0.375;
        case 5:    // 1/8(63x^5-70x^3 + 15*x )
            return 7.875 * x * x * x * x * x - 8.75 * x * x * x + 1.875 * x;
        case 6:    // 1/16(231x^6-315x^4+105x^2-5)
            return 14.4375 * x * x * x * x * x * x - 19.6875 * x * x * x * x + 6.5625 * x * x - 3.125;
        default: ErrorMessages::Error( "Legendre Polynomials only implemented up to order 6", CURRENT_FUNCTION ); return 0;
    }
}

void PNSolver::Solve() {
    auto log = spdlog::get( "event" );

    // angular flux at next time step (maybe store angular flux at all time steps, since time becomes energy?)
    VectorVector psiNew = _psi;
    double dFlux        = 1e10;
    Vector fluxNew( _nCells, 0.0 );
    Vector fluxOld( _nCells, 0.0 );

    for( unsigned i = 0; i < _nCells; ++i ) {
        _solverOutput[i] = _psi[i][0];
    }

    Save( -1 );    // Save initial condition

    VectorVector cellMids = _mesh->GetCellMidPoints();

    int rank;

    unsigned idx_system = 0;

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "dFlux" );

    // Loop over energies (pseudo-time of continuous slowing down approach)

    for( unsigned idx_energy = 0; idx_energy < _nEnergies; idx_energy++ ) {

        // Loop over all spatial cells
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;    // Dirichlet cells stay at IC, farfield assumption

            // Reset temporary variable psiNew
            for( int idx_lOrder = 0; idx_lOrder <= int( _nq ); idx_lOrder++ ) {
                for( int idx_kOrder = -idx_lOrder; idx_kOrder <= idx_lOrder; idx_kOrder++ ) {
                    idx_system                   = unsigned( GlobalIndex( idx_lOrder, idx_kOrder ) );
                    psiNew[idx_cell][idx_system] = 0.0;
                }
            }

            if( _mesh->GetDistanceToOrigin( idx_cell ) <= 0.1 ) {
                // std::cout << "Distance to origin: " << _mesh->GetDistanceToOrigin( idx_cell ) << "\n";
                // std::cout << "Old Cell value: \n " << _psi[idx_cell] << "\n" << psiNew[idx_cell] << "\n";
            }
            // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); idx_neighbor++ ) {

                // Compute flux contribution and store in psiNew to save memory
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                    _g->FluxVanLeer(
                        _Ax, _AxAbs, _Ay, _AyAbs, _Az, _AzAbs, _psi[idx_cell], _psi[idx_cell], _normals[idx_cell][idx_neighbor], psiNew[idx_cell] );
                // _g->Flux( _AxPlus,
                //           _AxMinus,
                //           _AyPlus,
                //           _AyMinus,
                //           _AzPlus,
                //           _AzMinus,
                //           _psi[idx_cell],
                //           _psi[idx_cell],
                //           _normals[idx_cell][idx_neighbor],
                //           psiNew[idx_cell] );
                else
                    _g->FluxVanLeer( _Ax,
                                     _AxAbs,
                                     _Ay,
                                     _AyAbs,
                                     _Az,
                                     _AzAbs,
                                     _psi[idx_cell],
                                     _psi[_neighbors[idx_cell][idx_neighbor]],
                                     _normals[idx_cell][idx_neighbor],
                                     psiNew[idx_cell] );
                //_g->Flux( _AxPlus,
                //          _AxMinus,
                //          _AyPlus,
                //          _AyMinus,
                //          _AzPlus,
                //          _AzMinus,
                //          _psi[idx_cell],
                //          _psi[_neighbors[idx_cell][idx_neighbor]],
                //          _normals[idx_cell][idx_neighbor],
                //          psiNew[idx_cell] );

                // if( _mesh->GetDistanceToOrigin( idx_cell ) <= 0.1 ) std::cout << "FluxPart " << idx_neighbor << " : \n" << psiNew[idx_cell] <<
                // "\n";
            }

            // time update angular flux with numerical flux and total scattering cross section
            for( int idx_lOrder = 0; idx_lOrder <= int( _nq ); idx_lOrder++ ) {
                for( int idx_kOrder = -idx_lOrder; idx_kOrder <= idx_lOrder; idx_kOrder++ ) {
                    idx_system = unsigned( GlobalIndex( idx_lOrder, idx_kOrder ) );

                    psiNew[idx_cell][idx_system] = _psi[idx_cell][idx_system] -
                                                   ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_system] /* cell averaged flux */
                                                   - _dE * _psi[idx_cell][idx_system] *
                                                         ( _sigmaA[idx_energy][idx_cell]                                    /* absorbtion influence */
                                                           + _sigmaS[idx_energy][idx_cell] * _scatterMatDiag[idx_system] ); /* scattering influence */
                }
            }
        }
        _psi = psiNew;

        double mass = 0.0;
        for( unsigned i = 0; i < _nCells; ++i ) {
            fluxNew[i]       = _psi[i][0];    // zeroth moment is raditation densitiy we are interested in
            _solverOutput[i] = _psi[i][0];
            mass += _psi[i][0];
        }

        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e} {:01.5e}", _energies[idx_energy], dFlux, mass );

        Save( idx_energy );
    }
}

void PNSolver::Save() const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<double> flux;
    flux.resize( _nCells );

    for( unsigned i = 0; i < _nCells; ++i ) {
        flux[i] = _psi[i][0];
    }
    std::vector<std::vector<double>> scalarField( 1, flux );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( _nEnergies ), results, fieldNames, _mesh );
}

void PNSolver::Save( int currEnergy ) const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<std::vector<double>> scalarField( 1, _solverOutput );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), results, fieldNames, _mesh );
}
