#include "solvers/pnsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/textprocessingtoolbox.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

PNSolver::PNSolver( Config* settings ) : Solver( settings ) {
    _LMaxDegree    = settings->GetMaxMomentDegree();
    _nTotalEntries = GlobalIndex( int( _LMaxDegree ), int( _LMaxDegree ) ) + 1;

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

    // std::cout << "System Matrix Set UP!" << std::endl;
    // Compute Decomposition in positive and negative (eigenvalue) parts of flux jacobians
    ComputeFluxComponents();

    // Compute diagonal of the scatter matrix (it's a diagonal matrix)
    ComputeScatterMatrix();

    // std::cout << "scatterMatrix : " << _scatterMatDiag << "\n";

    // AdaptTimeStep();

    if( settings->GetCleanFluxMat() ) CleanFluxMatrices();

    // std::cout << "--------\n";
    // std::cout << "_Ax :\n" << _Ax << "\n ";    // _AxP \n" << _AxPlus << "\n _AxM \n" << _AxMinus << "\n";
    // std::cout << "_Ay :\n" << _Ay << "\n ";    //_AyP \n" << _AyPlus << "\n _AyM \n" << _AyMinus << "\n";
    // std::cout << "_Az :\n" << _Az << "\n ";    //_AzP \n" << _AzPlus << "\n _AzM \n" << _AzMinus << "\n";
    //
    // std::cout << "_AxPlus :\n" << _AxPlus << "\n ";    // _AxP \n" << _AxPlus << "\n _AxM \n" << _AxMinus << "\n";
    // std::cout << "_AyPlus :\n" << _AyPlus << "\n ";    //_AyP \n" << _AyPlus << "\n _AyM \n" << _AyMinus << "\n";
    // std::cout << "_AzPlus :\n" << _AzPlus << "\n ";
    //
    // std::cout << "_AxMinus :\n" << _AxMinus << "\n ";    // _AxP \n" << _AxPlus << "\n _AxM \n" << _AxMinus << "\n";
    // std::cout << "_AyMinus :\n" << _AyMinus << "\n ";    //_AyP \n" << _AyPlus << "\n _AyM \n" << _AyMinus << "\n";
    // std::cout << "_AzMinus :\n" << _AzMinus << "\n ";
    //
    // std::cout << "_nCells: " << _nCells << "\n";

    // Solver output
    PrepareOutputFields();
}

void PNSolver::Solve() {

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    auto log = spdlog::get( "event" );

    VectorVector psiNew = _sol;

    double mass = 0;

    Save( -1 );    // Save initial condition

    unsigned idx_system = 0;

    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "mass" );

    // if( rank == 0 ) log->info( "{:03.8f}   {:01.5e} {:01.5e}", -1.0, dFlux, mass1 );

    // Loop over energies (pseudo-time of continuous slowing down approach)
    for( unsigned idx_energy = 0; idx_energy < _nEnergies; idx_energy++ ) {
        // Loop over all spatial cells
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;    // Dirichlet cells stay at IC, farfield assumption

            // Reset temporary variable psiNew
            for( int idx_lDegree = 0; idx_lDegree <= int( _LMaxDegree ); idx_lDegree++ ) {
                for( int idx_kOrder = -idx_lDegree; idx_kOrder <= idx_lDegree; idx_kOrder++ ) {
                    idx_system                   = unsigned( GlobalIndex( idx_lDegree, idx_kOrder ) );
                    psiNew[idx_cell][idx_system] = 0.0;
                }
            }

            // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
            for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); idx_neighbor++ ) {

                // Compute flux contribution and store in psiNew to save memory
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                    psiNew[idx_cell] += _g->Flux(
                        _AxPlus, _AxMinus, _AyPlus, _AyMinus, _AzPlus, _AzMinus, _sol[idx_cell], _sol[idx_cell], _normals[idx_cell][idx_neighbor] );
                else
                    psiNew[idx_cell] += _g->Flux( _AxPlus,
                                                  _AxMinus,
                                                  _AyPlus,
                                                  _AyMinus,
                                                  _AzPlus,
                                                  _AzMinus,
                                                  _sol[idx_cell],
                                                  _sol[_neighbors[idx_cell][idx_neighbor]],
                                                  _normals[idx_cell][idx_neighbor] );
            }

            // time update angular flux with numerical flux and total scattering cross section
            for( int idx_lOrder = 0; idx_lOrder <= int( _LMaxDegree ); idx_lOrder++ ) {
                for( int idx_kOrder = -idx_lOrder; idx_kOrder <= idx_lOrder; idx_kOrder++ ) {
                    idx_system = unsigned( GlobalIndex( idx_lOrder, idx_kOrder ) );

                    psiNew[idx_cell][idx_system] = _sol[idx_cell][idx_system] -
                                                   ( _dE / _areas[idx_cell] ) * psiNew[idx_cell][idx_system] /* cell averaged flux */
                                                   - _dE * _sol[idx_cell][idx_system] *
                                                         ( _sigmaA[idx_energy][idx_cell]                                    /* absorbtion influence */
                                                           + _sigmaS[idx_energy][idx_cell] * _scatterMatDiag[idx_system] ); /* scattering influence */
                }
            }
        }

        // Update Solution
        _sol = psiNew;

        // --- VTK and CSV Output ---
        mass = WriteOutputFields();
        Save( idx_energy );

        // --- Screen Output ---
        if( rank == 0 ) log->info( "{:03.8f}  {:01.5e}", _energies[idx_energy], mass );
    }
}

void PNSolver::ComputeSystemMatrices() {
    int idx_col      = 0;
    unsigned idx_row = 0;

    // loop over columns of A
    for( int idx_lOrder = 0; idx_lOrder <= int( _LMaxDegree ); idx_lOrder++ ) {          // index of legendre polynom
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
    // std::cout << "Eigenvalues x direction " << eigenValuesX << "\n";
    // std::cout << "Eigenvalues y direction " << eigenValuesY << "\n";
    // std::cout << "Eigenvalues z direction " << eigenValues << "\n";
    //
    // std::cout << "Spectral Radius X " << blaze::max( blaze::abs( eigenValuesX ) ) << "\n";
    // std::cout << "Spectral Radius Y " << blaze::max( blaze::abs( eigenValuesY ) ) << "\n";
    // std::cout << "Spectral Radius Z " << blaze::max( blaze::abs( eigenValues ) ) << "\n";

    // _combinedSpectralRadius = blaze::max( blaze::abs( eigenValues + eigenValuesX + eigenValuesY ) );
    // std::cout << "Spectral Radius combined " << _combinedSpectralRadius << "\n";
}

void PNSolver::ComputeScatterMatrix() {

    // --- Isotropic ---
    _scatterMatDiag[0] = 0.0;
    for( unsigned idx_diag = 1; idx_diag < _nTotalEntries; idx_diag++ ) {
        _scatterMatDiag[idx_diag] = 1.0;
    }
}

double PNSolver::LegendrePoly( double x, int l ) {    // Legacy. TO BE DELETED
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

void PNSolver::PrepareOutputFields() {
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

                _outputFields[0][0].resize( _nCells );
                _outputFieldNames[0][0] = "radiation flux density";
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
            default: ErrorMessages::Error( "Volume Output Group not defined for PN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

double PNSolver::WriteOutputFields() {
    double mass      = 0.0;
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();

    // Compute total "mass" of the system ==> to check conservation properties
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        mass += _sol[idx_cell][0] * _areas[idx_cell];    // Should probably go to postprocessing
    }

    for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
        switch( _settings->GetVolumeOutput()[idx_group] ) {
            case MINIMAL:
                for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                    _outputFields[idx_group][0][idx_cell] = _sol[idx_cell][0];
                }
                break;
            case MOMENTS:
                for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][idx_sys][idx_cell] = _sol[idx_cell][idx_sys];
                    }
                }
                break;
            default: ErrorMessages::Error( "Volume Output Group not defined for PN Solver!", CURRENT_FUNCTION ); break;
        }
    }
    return mass;
}

void PNSolver::Save() const { ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh ); }

void PNSolver::Save( int currEnergy ) const {
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), _outputFields, _outputFieldNames, _mesh );
}

void PNSolver::CleanFluxMatrices() {
    for( unsigned idx_row = 0; idx_row < _nTotalEntries; idx_row++ ) {
        for( unsigned idx_col = 0; idx_col < _nTotalEntries; idx_col++ ) {
            if( std::abs( _AxAbs( idx_row, idx_col ) ) < 0.00000000001 ) _AxAbs( idx_row, idx_col ) = 0.0;
            if( std::abs( _AxPlus( idx_row, idx_col ) ) < 0.00000000001 ) _AxPlus( idx_row, idx_col ) = 0.0;
            if( std::abs( _AxMinus( idx_row, idx_col ) ) < 0.00000000001 ) _AxMinus( idx_row, idx_col ) = 0.0;

            if( std::abs( _AyAbs( idx_row, idx_col ) ) < 0.00000000001 ) _AyAbs( idx_row, idx_col ) = 0.0;
            if( std::abs( _AyPlus( idx_row, idx_col ) ) < 0.00000000001 ) _AyPlus( idx_row, idx_col ) = 0.0;
            if( std::abs( _AyMinus( idx_row, idx_col ) ) < 0.00000000001 ) _AyMinus( idx_row, idx_col ) = 0.0;

            if( std::abs( _AzAbs( idx_row, idx_col ) ) < 0.00000000001 ) _AzAbs( idx_row, idx_col ) = 0.0;
            if( std::abs( _AzPlus( idx_row, idx_col ) ) < 0.00000000001 ) _AzPlus( idx_row, idx_col ) = 0.0;
            if( std::abs( _AzMinus( idx_row, idx_col ) ) < 0.00000000001 ) _AzMinus( idx_row, idx_col ) = 0.0;
        }
    }
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
    if( l >= 0 && l <= int( _LMaxDegree ) ) {
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
