#include "solvers/pnsolver.h"
#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "fluxes/numericalflux.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/textprocessingtoolbox.h"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

PNSolver::PNSolver( Config* settings ) : Solver( settings ) {
    _LMaxDegree    = settings->GetMaxMomentDegree();
    _nTotalEntries = GlobalIndex( int( _LMaxDegree ), int( _LMaxDegree ) ) + 1;

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

    _solDx = VectorVector( _nCells, Vector( _nTotalEntries, 0.0 ) );
    _solDy = VectorVector( _nCells, Vector( _nTotalEntries, 0.0 ) );

    // Fill System Matrices
    ComputeSystemMatrices();

    // Compute Decomposition in positive and negative (eigenvalue) parts of flux jacobians
    ComputeFluxComponents();

    // Compute diagonal of the scatter matrix (it's a diagonal matrix)
    ComputeScatterMatrix();

    // AdaptTimeStep();

    if( settings->GetCleanFluxMat() ) CleanFluxMatrices();

    // Compute moments of initial condition
    // TODO

    // Solver output
    PrepareVolumeOutput();
}

void PNSolver::IterPreprocessing() {
    // Nothing to preprocess for PNSolver
}

void PNSolver::IterPostprocessing() {
    // --- Update Solution ---
    _sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();
}

void PNSolver::ComputeRadFlux() {
    double firstMomentScaleFactor = sqrt( 4 * M_PI );
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        _fluxNew[idx_cell] = _sol[idx_cell][0] * firstMomentScaleFactor;
    }
}

void PNSolver::FluxUpdate() {
    if( _reconsOrder > 1 ) {
        _mesh->ReconstructSlopesU( _nTotalEntries, _solDx, _solDy, _sol );    // unstructured reconstruction
        //_mesh->ComputeSlopes( _nTotalEntries, _solDx, _solDy, _sol );    // unstructured reconstruction
    }
    //Vector solL( _nTotalEntries );
    //Vector solR( _nTotalEntries );
    auto solL = _sol[2];
    auto solR = _sol[2];

    // Loop over all spatial cells
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {

        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;

        // Reset temporary variable psiNew
        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            _solNew[idx_cell][idx_sys] = 0.0;
        }

        // Loop over all neighbor cells (edges) of cell j and compute numerical fluxes
        for( unsigned idx_neighbor = 0; idx_neighbor < _neighbors[idx_cell].size(); idx_neighbor++ ) {

            // Compute flux contribution and store in psiNew to save memory
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_neighbor] == _nCells )
                _solNew[idx_cell] += _g->Flux(
                    _AxPlus, _AxMinus, _AyPlus, _AyMinus, _AzPlus, _AzMinus, _sol[idx_cell], _sol[idx_cell], _normals[idx_cell][idx_neighbor] );
            else {
                switch( _reconsOrder ) {
                    // first order solver
                    case 1:
                        _solNew[idx_cell] += _g->Flux( _AxPlus,
                                                       _AxMinus,
                                                       _AyPlus,
                                                       _AyMinus,
                                                       _AzPlus,
                                                       _AzMinus,
                                                       _sol[idx_cell],
                                                       _sol[_neighbors[idx_cell][idx_neighbor]],
                                                       _normals[idx_cell][idx_neighbor] );
                        break;
                    // second order solver
                    case 2:
                        // left status of interface
                        solL = _sol[idx_cell] +
                               _solDx[idx_cell] * ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[idx_cell][0] ) + 
                               _solDy[idx_cell] * ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[idx_cell][1] );
                        // right status of interface
                        solR = _sol[_neighbors[idx_cell][idx_neighbor]] +
                               _solDx[_neighbors[idx_cell][idx_neighbor]] * ( _interfaceMidPoints[idx_cell][idx_neighbor][0] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][0] ) + 
                               _solDy[_neighbors[idx_cell][idx_neighbor]] * ( _interfaceMidPoints[idx_cell][idx_neighbor][1] - _cellMidPoints[_neighbors[idx_cell][idx_neighbor]][1] );
                        // positivity checker (if not satisfied, deduce to first order)
                        // I manually turned it off here since Pn produces negative solutions essentially
                        //if( min(solL) < 0.0 || min(solR) < 0.0 ) {
                        //    solL = _sol[idx_cell];
                        //    solR = _sol[_neighbors[idx_cell][idx_neighbor]];
                        //}
                        // flux evaluation
                        _solNew[idx_cell] += _g->Flux( _AxPlus,
                                                       _AxMinus,
                                                       _AyPlus,
                                                       _AyMinus,
                                                       _AzPlus,
                                                       _AzMinus,
                                                       solL,
                                                       solR,
                                                       _normals[idx_cell][idx_neighbor] );
                        break;
                    // default: first order solver
                    default: 
                        _solNew[idx_cell] += _g->Flux( _AxPlus,
                                                       _AxMinus,
                                                       _AyPlus,
                                                       _AyMinus,
                                                       _AzPlus,
                                                       _AzMinus,
                                                       _sol[idx_cell],
                                                       _sol[_neighbors[idx_cell][idx_neighbor]],
                                                       _normals[idx_cell][idx_neighbor] );
                }
            }
        }
    }
}

void PNSolver::FVMUpdate( unsigned idx_energy ) {
    // Loop over all spatial cells
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        // Dirichlet cells stay at IC, farfield assumption
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // Flux update
        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            _solNew[idx_cell][idx_sys] = _sol[idx_cell][idx_sys] - ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_sys] /* cell averaged flux */
                                         - _dE * _sol[idx_cell][idx_sys] *
                                               ( _sigmaT[idx_energy][idx_cell]                                 /* absorbtion influence */
                                                 + _sigmaS[idx_energy][idx_cell] * _scatterMatDiag[idx_sys] ); /* scattering influence */
        }
        // Source Term
        _solNew[idx_cell][0] += _dE * _Q[0][idx_cell][0];
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

        // eigenValuesX = eigenValues;
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

        // eigenValuesY = eigenValues;
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
    _scatterMatDiag[0] = -1.0;
    for( unsigned idx_diag = 1; idx_diag < _nTotalEntries; idx_diag++ ) {
        _scatterMatDiag[idx_diag] = 0.0;
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

void PNSolver::PrepareVolumeOutput() {
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
            default: ErrorMessages::Error( "Volume Output Group not defined for PN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

void PNSolver::WriteVolumeOutput( unsigned idx_pseudoTime ) {
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
    }
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
