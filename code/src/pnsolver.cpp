
#include "pnsolver.h"
#include "../ext/blaze/blaze/math/blas/cblas/gemm.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/textprocessingtoolbox.h"
#include <mpi.h>

PNSolver::PNSolver( Config* settings ) : Solver( settings ) {
    _nTotalEntries = GlobalIndex( _nq, _nq ) + 1;

    // transform sigmaT and sigmaS in sigmaA.
    _sigmaA = VectorVector( _nEnergies, Vector( _nCells ) );    // Get rid of this extra vektor!
    for( unsigned n = 0; n < _nEnergies; n++ ) {
        for( unsigned j = 0; j < _nCells; j++ ) {
            _sigmaA[n][j] = _sigmaT[n][j] - _sigmaS[n][j];
        }
    }

    // Initialize System Matrices
    _Ax = SymMatrix( _nTotalEntries );
    _Ay = SymMatrix( _nTotalEntries );
    _Az = SymMatrix( _nTotalEntries );

    _AxPlus  = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AxMinus = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AyPlus  = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AyMinus = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AzPlus  = Matrix( _nTotalEntries, _nTotalEntries, 0 );
    _AzMinus = Matrix( _nTotalEntries, _nTotalEntries, 0 );

    // Fill System Matrices
    SetupSystemMatrices();

    std::cout << "System Matrix Set UP!" << std::endl;

    // Compute Eigenvalues of the system matrices
    ComputeFluxComponents();
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

unsigned PNSolver::GlobalIndex( int l, int k ) const {
    int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
    int prevIndicesThisLevel = k + l;    // number of previous indices in current level
    return unsigned( numIndicesPrevLevel + prevIndicesThisLevel );
}

int PNSolver::Sgn( int k ) const {
    if( k >= 0 )
        return 1;
    else
        return -1;
}

void PNSolver::SetupSystemMatrices() {
    int j      = 0;
    unsigned i = 0;

    // loop over columns of A
    for( int l = 0; l <= int( _nq ); ++l ) {    // index of legendre polynom
        for( int k = -l; k <= l; ++k ) {        // second index of legendre function
            i = unsigned( GlobalIndex( l, k ) );

            // flux matrix in direction x
            if( k != -1 ) {
                j = GlobalIndex( l - 1, kMinus( k ) );
                if( j >= 0 && j < int( _nTotalEntries ) ) {
                    _Ax( i, unsigned( j ) ) = 0.5 * CTilde( l - 1, std::abs( k ) - 1 );
                }

                j = GlobalIndex( l + 1, kMinus( k ) );
                if( j >= 0 && j < int( _nTotalEntries ) ) {
                    _Ax( i, unsigned( j ) ) = -0.5 * DTilde( l + 1, std::abs( k ) - 1 );
                }
            }

            j = GlobalIndex( l - 1, kPlus( k ) );
            if( j >= 0 && j < int( _nTotalEntries ) ) {
                _Ax( i, unsigned( j ) ) = -0.5 * ETilde( l - 1, std::abs( k ) + 1 );
            }

            j = GlobalIndex( l + 1, kPlus( k ) );
            if( j >= 0 && j < int( _nTotalEntries ) ) {
                _Ax( i, unsigned( j ) ) = 0.5 * FTilde( l + 1, std::abs( k ) + 1 );
            }

            //
            // flux matrix in direction y
            if( k != 1 ) {
                j = GlobalIndex( l - 1, -kMinus( k ) );
                if( j >= 0 && j < int( _nTotalEntries ) ) {
                    _Ay( i, unsigned( j ) ) = -0.5 * Sgn( k ) * CTilde( l - 1, std::abs( k ) - 1 );
                }

                j = GlobalIndex( l + 1, -kMinus( k ) );
                if( j >= 0 && j < int( _nTotalEntries ) ) {
                    _Ay( i, unsigned( j ) ) = 0.5 * Sgn( k ) * DTilde( l + 1, std::abs( k ) - 1 );
                }
            }

            j = GlobalIndex( l - 1, -kPlus( k ) );
            if( j >= 0 && j < int( _nTotalEntries ) ) {
                _Ay( i, unsigned( j ) ) = -0.5 * Sgn( k ) * ETilde( l - 1, std::abs( k ) + 1 );
            }

            j = GlobalIndex( l + 1, -kPlus( k ) );
            if( j >= 0 && j < int( _nTotalEntries ) ) {
                _Ay( i, unsigned( j ) ) = 0.5 * Sgn( k ) * FTilde( l + 1, std::abs( k ) + 1 );
            }

            //
            // flux matrix in direction z
            j = GlobalIndex( l - 1, k );
            if( j >= 0 && j < int( _nTotalEntries ) ) {
                _Az( i, unsigned( j ) ) = AParam( l - 1, k );
            }

            j = GlobalIndex( l + 1, k );
            if( j >= 0 && j < int( _nTotalEntries ) ) {
                _Az( i, unsigned( j ) ) = BParam( l + 1, k );
            }
        }
    }
}

void PNSolver::ComputeFluxComponents() {
    Vector eigenValues( _nTotalEntries );
    MatrixCol eigenVectors( _nTotalEntries, _nTotalEntries, 0 );    // ColumnMatrix for _AxPlus * eigenVectors Multiplication via SIMD
    // --- For x Direction ---
    {
        std::cout << "A_x = \n" << _Ax << std::endl;

        blaze::eigen( _Ax, eigenValues, eigenVectors );    // Compute Eigenvalues and Eigenvectors

        // Compute Flux Matrices A+ and A-
        for( unsigned idx_ij = 0; idx_ij < _nTotalEntries; idx_ij++ ) {
            if( eigenValues[idx_ij] >= 0 ) {
                _AxPlus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // positive part of Diagonal Matrix stored in _AxPlus
            }
            else {
                _AxMinus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // negative part of Diagonal Matrix stored in _AxMinus
            }
        }

        _AxPlus  = eigenVectors * _AxPlus;    // col*row minimum performance
        _AxMinus = eigenVectors * _AxMinus;
        blaze::transpose( eigenVectors );
        _AxPlus  = _AxPlus * eigenVectors;    // row*col maximum performance
        _AxMinus = _AxMinus * eigenVectors;

        std::cout << " _AxPlus = \n" << _AxPlus << "\n";
        std::cout << " _AxMinus = \n" << _AxMinus << "\n";
        std::cout << "A_x = \n" << _AxPlus + _AxMinus << std::endl;
    }
    // --- For y Direction -------
    {
        blaze::eigen( _Ay, eigenValues, eigenVectors );    // Compute Eigenvalues and Eigenvectors

        // Compute Flux Matrices A+ and A-
        for( unsigned idx_ij = 0; idx_ij < _nTotalEntries; idx_ij++ ) {
            if( eigenValues[idx_ij] >= 0 ) {
                _AyPlus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // positive part of Diagonal Matrix stored in _AxPlus
            }
            else {
                _AyMinus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // negative part of Diagonal Matrix stored in _AxMinus
            }
        }

        _AyPlus  = eigenVectors * _AyPlus;
        _AyMinus = eigenVectors * _AyMinus;
        blaze::transpose( eigenVectors );
        _AyPlus  = _AyPlus * eigenVectors;
        _AyMinus = _AyMinus * eigenVectors;
    }
    // --- For z Direction -------
    {
        blaze::eigen( _Az, eigenValues, eigenVectors );    // Compute Eigenvalues and Eigenvectors

        // Compute Flux Matrices A+ and A-
        for( unsigned idx_ij = 0; idx_ij < _nTotalEntries; idx_ij++ ) {
            if( eigenValues[idx_ij] >= 0 ) {
                _AzPlus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // positive part of Diagonal Matrix stored in _AxPlus
            }
            else {
                _AzMinus( idx_ij, idx_ij ) = eigenValues[idx_ij];    // negative part of Diagonal Matrix stored in _AxMinus
            }
        }

        _AzPlus  = eigenVectors * _AzPlus;
        _AzMinus = eigenVectors * _AzMinus;
        blaze::transpose( eigenVectors );
        _AzPlus  = _AzPlus * eigenVectors;
        _AzMinus = _AzMinus * eigenVectors;
    }
}

double PNSolver::ComputeScatterMatrix( int l ) {
    // right now independent of spatial coordinates
    // G[i][i] = 1- g_l
    // g_l = 2*pi*\int_-1 ^1 k(x,my)P_l(my)dmy
    // use midpoint rule.
    // Here is room for optimization!

    // Only isotropic scattering: k = 1/(2pi) in 2d
    unsigned nSteps    = 50;
    double integralSum = 0.5 * ( Legendre( -1.0, l ) + Legendre( 1.0, l ) );

    for( unsigned i = 1; i < nSteps - 1; i++ ) {
        integralSum += Legendre( -1.0 + double( i ) * 2 / double( nSteps ), l );
    }
    integralSum *= 2 / double( nSteps );

    // isotropic scattering and prefactor of eigenvalue
    integralSum *= 0.5;    // 2pi/4pi

    return 1 - integralSum;
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

    int rank;

    int idx_system = 0;

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "dFlux" );

    // loop over energies (pseudo-time)
    for( unsigned idx_energy = 0; idx_energy < _nEnergies; ++idx_energy ) {
        // loop over all spatial cells
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
            // loop over all ordinates
            for( unsigned idx_lOrder = 0; idx_lOrder < _nq; ++idx_lOrder ) {
                for( int idx_kOrder = -idx_lOrder; idx_kOrder < int( idx_lOrder ); ++idx_lOrder ) {

                    idx_system = GlobalIndex( idx_lOrder, idx_kOrder );

                    psiNew[idx_cell][idx_system] = 0.0;

                    // loop over all neighbor cells (edges) of cell j and compute numerical fluxes
                    for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); ++idx_neighbor ) {

                        // store flux contribution on psiNew to save memory
                        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                            psiNew[idx_cell][idx_system] += _g->Flux( _quadPoints[idx_system],
                                                                      _psi[idx_cell][idx_system],
                                                                      _psi[idx_cell][idx_system],
                                                                      _normals[idx_cell][idx_neighbor] );    // TODO Flux computation!
                        else
                            psiNew[idx_cell][idx_system] += _g->Flux( _quadPoints[idx_system],
                                                                      _psi[idx_cell][idx_system],
                                                                      _psi[_neighbors[idx_cell][idx_neighbor]][idx_system],
                                                                      _normals[idx_cell][idx_neighbor] );    // TODO Flux computation!
                    }
                    // time update angular flux with numerical flux and total scattering cross section
                    psiNew[idx_cell][idx_system] =
                        _psi[idx_cell][idx_system] - ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_system] -
                        _dE * ( _sigmaA[idx_energy][idx_cell] + _sigmaS[idx_energy][idx_cell] * ComputeScatterMatrix( int( idx_system ) ) ) *
                            _psi[idx_cell][idx_system];
                }
            }
            // compute scattering effects
            // psiNew[idx_cell] += _dE * _sigmaS[idx_energy][idx_cell] * _scatteringKernel * _psi[idx_cell];    // multiply scattering matrix with psi

            // TODO: Adapt to PN
            // TODO: figure out a more elegant way
            // add external source contribution
            if( _Q.size() == 1u ) {                   // constant source for all energies
                if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                    psiNew[idx_cell] += _dE * _Q[0][idx_cell][0];
                else
                    psiNew[idx_cell] += _dE * _Q[0][idx_cell];
            }
            else {
                if( _Q[0][idx_cell].size() == 1u )    // isotropic source
                    psiNew[idx_cell] += _dE * _Q[idx_energy][idx_cell][0];
                else
                    psiNew[idx_cell] += _dE * _Q[idx_energy][idx_cell];
            }
        }
        _psi = psiNew;
        for( unsigned i = 0; i < _nCells; ++i ) {
            fluxNew[i] = dot( _psi[i], _weights );
        }
        dFlux   = blaze::l2Norm( fluxNew - fluxOld );
        fluxOld = fluxNew;
        if( rank == 0 ) log->info( "{:03.8f}   {:01.5e}", _energies[idx_energy], dFlux );
    }
}

void PNSolver::Save() const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<double> flux( _nCells, 0.0 );
    for( unsigned i = 0; i < _nCells; ++i ) {
        flux[i] = dot( _psi[i], _weights );
    }
    std::vector<std::vector<double>> scalarField( 1, flux );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile(), results, fieldNames, _mesh );
}
