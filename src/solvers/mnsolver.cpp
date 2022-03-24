#include "solvers/mnsolver.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "entropies/entropybase.hpp"
#include "fluxes/numericalflux.hpp"
#include "optimizers/optimizerbase.hpp"
#include "problems/problembase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "velocitybasis/sphericalbase.hpp"

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

MNSolver::MNSolver( Config* settings ) : SolverBase( settings ) {

    _polyDegreeBasis = settings->GetMaxMomentDegree();
    _basis           = SphericalBase::Create( _settings );
    _nSystem         = _basis->GetBasisSize();

    // build quadrature object and store quadrature points and weights
    _quadPoints = _quadrature->GetPoints();
    _weights    = _quadrature->GetWeights();
    //_nq               = _quadrature->GetNq();
    _quadPointsSphere = _quadrature->GetPointsSphere();
    //_settings->SetNQuadPoints( _nq );

    // Initialize Entropy
    _entropy = EntropyBase::Create( _settings );

    // Initialize Optimizer
    _optimizer = OptimizerBase::Create( _settings );

    // Initialize lagrange Multiplier
    _alpha = VectorVector( _nCells, Vector( _nSystem, 0.0 ) );

    // Initialize kinetic density at grid cells
    _kineticDensity = VectorVector( _nCells, Vector( _nq, 0.0 ) );

    // Limiter variables
    _solDx   = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    _solDy   = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    _limiter = VectorVector( _nCells, Vector( _nq, 0.0 ) );

    // Initialize and Pre-Compute Moments at quadrature points
    _momentBasis = VectorVector( _nq, Vector( _nSystem, 0.0 ) );
    ComputeMoments();

    // Initialize Scatter Matrix --
    _scatterMatDiag = Vector( _nSystem, 0.0 );
    ComputeScatterMatrix();
}

MNSolver::~MNSolver() {
    delete _entropy;
    delete _optimizer;
    delete _basis;
}

void MNSolver::ComputeScatterMatrix() {

    // --- Isotropic ---
    if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
        _scatterMatDiag[0] = 1.0;
        for( unsigned idx_diag = 1; idx_diag < _nSystem; idx_diag++ ) {
            _scatterMatDiag[idx_diag] = 0.0;    // SPHERICAL_HARMONICS are orthogonal basis
        }
    }
    else {    // SPHERICAL_MONOMIALS
        for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
            _scatterMatDiag[idx_sys] = 0.0;
            for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                _scatterMatDiag[idx_sys] += _momentBasis[idx_quad][idx_sys] * _weights[idx_quad];
            }
            if( _settings->GetDim() == 1 ) {
                _scatterMatDiag[idx_sys] /= 2;
            }
            else if( _settings->GetDim() == 2 ) {
                _scatterMatDiag[idx_sys] /= M_PI;
            }
            else {    // 3D
                _scatterMatDiag[idx_sys] /= 4.0 * M_PI;
            }
        }
    }
}

void MNSolver::ComputeMoments() {
    double my, phi;

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my                     = _quadPointsSphere[idx_quad][0];
        phi                    = _quadPointsSphere[idx_quad][1];
        _momentBasis[idx_quad] = _basis->ComputeSphericalBasis( my, phi );
    }
}

Vector MNSolver::ConstructFlux( unsigned idx_cell ) {
    //--- Integration of moments of flux ---
    // Delete here and in csdmn solver
}

void MNSolver::ComputeRealizableSolution( unsigned idx_cell ) {
    // double entropyReconstruction = 0.0;

    for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {    // reset solution
        _sol[idx_cell][idx_sys] = 0.0;
    }
    _optimizer->ReconstructMoments( _sol[idx_cell], _alpha[idx_cell], _momentBasis );
}

void MNSolver::IterPreprocessing( unsigned /*idx_pseudotime*/ ) {

    // ------- Entropy closure Step ----------------
    _optimizer->SolveMultiCell( _alpha, _sol, _momentBasis );    // parallel

    // ------- Solution reconstruction step ----
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            // compute the kinetic density at all grid cells
            _kineticDensity[idx_cell][idx_quad] = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_cell], _momentBasis[idx_quad] ) );
        }
        if( _settings->GetRealizabilityReconstruction() ) ComputeRealizableSolution( idx_cell );
    }
    // ------ Compute slope limiters and cell gradients ---
    if( _reconsOrder > 1 ) {
        _mesh->ComputeSlopes( _nq, _solDx, _solDy, _kineticDensity );               // parallel
        _mesh->ComputeLimiter( _nq, _solDx, _solDy, _kineticDensity, _limiter );    // parallel
    }
}

void MNSolver::IterPostprocessing( unsigned /*idx_iter*/ ) {
    // --- Update Solution ---
    //_sol = _solNew;

    // --- Compute Flux for solution and Screen Output ---
    ComputeRadFlux();
}

void MNSolver::ComputeRadFlux() {
    double firstMomentScaleFactor = 4 * M_PI;
    if( _settings->GetProblemName() == PROBLEM_Aircavity1D || _settings->GetProblemName() == PROBLEM_Linesource1D ||
        _settings->GetProblemName() == PROBLEM_Checkerboard1D ) {
        firstMomentScaleFactor = 2.0;
    }
    if( _settings->GetDim() == 2 ) {
        firstMomentScaleFactor = M_PI;
    }
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        _fluxNew[idx_cell] = _sol[idx_cell][0] * firstMomentScaleFactor;
    }
}

void MNSolver::FluxUpdate() {
    // Loop over all spatial cells
    if( _settings->GetProblemName() == PROBLEM_Aircavity1D || _settings->GetProblemName() == PROBLEM_Linesource1D ||
        _settings->GetProblemName() == PROBLEM_Checkerboard1D ) {
        FluxUpdatePseudo1D();
    }
    else {
        FluxUpdatePseudo2D();
    }
}

void MNSolver::FluxUpdatePseudo1D() {
// Loop over the grid cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        // Dirichlet Boundaries stayd
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        double solL, solR, kineticFlux;
        Vector flux( _nSystem, 0.0 );
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            kineticFlux = 0.0;    // reset temorary flux

            for( unsigned idx_nbr = 0; idx_nbr < _neighbors[idx_cell].size(); idx_nbr++ ) {
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_nbr] == _nCells ) {
                    // Boundary cells are first order and mirror ghost cells
                    solL = _kineticDensity[idx_cell][idx_quad];
                    solR = solL;
                }
                else {
                    // interior cell
                    unsigned int nbr_glob = _neighbors[idx_cell][idx_nbr];    // global idx of neighbor cell
                    if( _reconsOrder == 1 ) {
                        solL = _kineticDensity[idx_cell][idx_quad];
                        solR = _kineticDensity[_neighbors[idx_cell][idx_nbr]][idx_quad];
                    }
                    else if( _reconsOrder == 2 ) {
                        solL = _kineticDensity[idx_cell][idx_quad] +
                               _limiter[idx_cell][idx_quad] *
                                   ( _solDx[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[idx_cell][0] ) +
                                     _solDy[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[idx_cell][1] ) );
                        solR = _kineticDensity[nbr_glob][idx_quad] +
                               _limiter[nbr_glob][idx_quad] *
                                   ( _solDx[nbr_glob][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[nbr_glob][0] ) +
                                     _solDy[nbr_glob][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[nbr_glob][1] ) );
                    }
                    else {
                        ErrorMessages::Error( "Reconstruction order not supported.", CURRENT_FUNCTION );
                    }
                }
                // Kinetic flux
                kineticFlux += _g->Flux1D( _quadPoints[idx_quad], solL, solR, _normals[idx_cell][idx_nbr] );
            }
            // Solution flux
            flux += _momentBasis[idx_quad] * ( _weights[idx_quad] * kineticFlux );
        }
        _solNew[idx_cell] = flux;
    }
}

void MNSolver::FluxUpdatePseudo2D() {
// Loop over the grid cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        // Dirichlet Boundaries stayd
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        double solL, solR, kineticFlux;
        Vector flux( _nSystem, 0.0 );
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            kineticFlux = 0.0;    // reset temorary flux
            for( unsigned idx_nbr = 0; idx_nbr < _neighbors[idx_cell].size(); idx_nbr++ ) {
                if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[idx_cell][idx_nbr] == _nCells ) {
                    // Boundary cells are first order and mirror ghost cells
                    solL = _kineticDensity[idx_cell][idx_quad];
                    solR = solL;
                }
                else {
                    // interior cell
                    unsigned int nbr_glob = _neighbors[idx_cell][idx_nbr];    // global idx of neighbor cell
                    if( _reconsOrder == 1 ) {
                        solL = _kineticDensity[idx_cell][idx_quad];
                        solR = _kineticDensity[_neighbors[idx_cell][idx_nbr]][idx_quad];
                    }
                    else if( _reconsOrder == 2 ) {
                        solL = _kineticDensity[idx_cell][idx_quad] +
                               _limiter[idx_cell][idx_quad] *
                                   ( _solDx[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[idx_cell][0] ) +
                                     _solDy[idx_cell][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[idx_cell][1] ) );
                        solR = _kineticDensity[nbr_glob][idx_quad] +
                               _limiter[nbr_glob][idx_quad] *
                                   ( _solDx[nbr_glob][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[nbr_glob][0] ) +
                                     _solDy[nbr_glob][idx_quad] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[nbr_glob][1] ) );
                    }
                    else {
                        ErrorMessages::Error( "Reconstruction order not supported.", CURRENT_FUNCTION );
                    }
                }
                // Kinetic flux
                kineticFlux += _g->Flux( _quadPoints[idx_quad], solL, solR, _normals[idx_cell][idx_nbr] );
            }
            // Solution flux
            flux += _momentBasis[idx_quad] * ( _weights[idx_quad] * kineticFlux );
        }
        _solNew[idx_cell] = flux;
    }
}

void MNSolver::FVMUpdate( unsigned idx_iter ) {
// Loop over the grid cells
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        // Dirichlet Boundaries stay
        if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) continue;
        // Flux update
        for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
            _solNew[idx_cell][idx_sys] = _sol[idx_cell][idx_sys]                                                            /* old solution */
                                         - ( _dE / _areas[idx_cell] ) * _solNew[idx_cell][idx_sys]                          /* cell averaged flux */
                                         + _dE * _sigmaS[idx_iter][idx_cell] * _sol[idx_cell][0] * _scatterMatDiag[idx_sys] /* scattering gain */
                                         - _dE * ( _sigmaT[idx_iter][idx_cell] ) * _sol[idx_cell][idx_sys]; /* scattering and absorbtion loss */
        }
        // Source Term
        _solNew[idx_cell] += _dE * _Q[0][idx_cell];
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
                _outputFields[idx_group].resize( _nSystem );
                _outputFieldNames[idx_group].resize( _nSystem );
                if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
                    if( _settings->GetDim() == 3 ) {
                        for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {    // 3D
                            for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                                _outputFields[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )].resize( _nCells );
                                _outputFieldNames[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )] =
                                    std::string( "u_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                            }
                        }
                    }
                    else if( _settings->GetDim() == 2 ) {
                        unsigned count = 0;
                        for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                            for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                                if( idx_l == 0 || idx_l != idx_k ) {
                                    _outputFields[idx_group][count].resize( _nCells );
                                    _outputFieldNames[idx_group][count] =
                                        std::string( "u_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                                    count++;
                                }
                            }
                        }
                    }
                    else if( _settings->GetDim() == 1 ) {
                        for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                            _outputFields[idx_group][idx_l].resize( _nCells );
                            _outputFieldNames[idx_group][idx_l] = std::string( "u_" + std::to_string( idx_l ) + "^0" );
                        }
                    }
                }
                else {
                    for( unsigned idx_l = 0; idx_l <= _polyDegreeBasis; idx_l++ ) {
                        unsigned maxOrder_k = _basis->GetCurrDegreeSize( idx_l );
                        for( unsigned idx_k = 0; idx_k < maxOrder_k; idx_k++ ) {
                            _outputFields[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )].resize( _nCells );
                            _outputFieldNames[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )] =
                                std::string( "u_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                        }
                    }
                }
                break;
            case DUAL_MOMENTS:
                // As many entries as there are moments in the system
                _outputFields[idx_group].resize( _nSystem );
                _outputFieldNames[idx_group].resize( _nSystem );
                if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
                    if( _settings->GetDim() == 3 ) {
                        for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                            for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                                _outputFields[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )].resize( _nCells );
                                _outputFieldNames[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )] =
                                    std::string( "alpha_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                            }
                        }
                    }
                    else if( _settings->GetDim() == 2 ) {
                        unsigned count = 0;
                        for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                            for( int idx_k = -idx_l; idx_k <= idx_l; idx_k++ ) {
                                if( idx_l == 0 || idx_l != idx_k ) {
                                    _outputFields[idx_group][count].resize( _nCells );
                                    _outputFieldNames[idx_group][count] =
                                        std::string( "alpha_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                                    count++;
                                }
                            }
                        }
                    }
                    else if( _settings->GetDim() == 1 ) {
                        for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                            _outputFields[idx_group][idx_l].resize( _nCells );
                            _outputFieldNames[idx_group][idx_l] = std::string( "alpha_" + std::to_string( idx_l ) + "^0" );
                        }
                    }
                }
                else {    // SPHERICAL_MONOMIALS
                    for( int idx_l = 0; idx_l <= (int)_polyDegreeBasis; idx_l++ ) {
                        unsigned maxOrder_k = _basis->GetCurrDegreeSize( idx_l );
                        for( unsigned idx_k = 0; idx_k < maxOrder_k; idx_k++ ) {
                            _outputFields[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )].resize( _nCells );
                            _outputFieldNames[idx_group][_basis->GetGlobalIndexBasis( idx_l, idx_k )] =
                                std::string( "alpha_" + std::to_string( idx_l ) + "^" + std::to_string( idx_k ) );
                        }
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

void MNSolver::WriteVolumeOutput( unsigned idx_iter ) {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();
    // Check if volume output fields are written to file this iteration
    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_iter == _nEnergies - 1 ) /* need sol at last iteration */ ) {
        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _fluxNew[idx_cell];
                    }
                    break;
                case MOMENTS:
                    for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
                        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                            _outputFields[idx_group][idx_sys][idx_cell] = _sol[idx_cell][idx_sys];
                        }
                    }
                    break;
                case DUAL_MOMENTS:
                    for( unsigned idx_sys = 0; idx_sys < _nSystem; idx_sys++ ) {
                        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
                            _outputFields[idx_group][idx_sys][idx_cell] = _alpha[idx_cell][idx_sys];
                        }
                    }
                    break;
                case ANALYTIC:
                    // Compute total "mass" of the system ==> to check conservation properties
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        double time                           = idx_iter * _dE;
                        _outputFields[idx_group][0][idx_cell] = _problem->GetAnalyticalSolution(
                            _mesh->GetCellMidPoints()[idx_cell][0], _mesh->GetCellMidPoints()[idx_cell][1], time, _sigmaS[idx_iter][idx_cell] );
                    }
                    break;
                default: ErrorMessages::Error( "Volume Output Group not defined for MN Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}
