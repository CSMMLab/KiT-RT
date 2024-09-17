#ifdef BUILD_MPI
#include <mpi.h>
#endif
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "kernels/scatteringkernelbase.hpp"
#include "problems/problembase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "solvers/snsolver_hpc.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

SNSolverHPC::SNSolverHPC( Config* settings ) {
#ifdef BUILD_MPI
    // Initialize MPI
    MPI_Comm_size( MPI_COMM_WORLD, &_numProcs );
    MPI_Comm_rank( MPI_COMM_WORLD, &_rank );
#endif
#ifndef BUILD_MPI
    _numProcs = 1;    // default values
    _rank     = 0;
#endif

    _settings = settings;
    _currTime = 0.0;

    _nOutputMoments = 2;    //  Currently only u_1 (x direction) and u_1 (y direction) are supported
    // Create Mesh
    _mesh = LoadSU2MeshFromFile( settings );
    _settings->SetNCells( _mesh->GetNumCells() );
    auto quad = QuadratureBase::Create( settings );
    _settings->SetNQuadPoints( quad->GetNq() );

    _nCells = _mesh->GetNumCells();
    _nNbr   = _mesh->GetNumNodesPerCell();
    _nDim   = _mesh->GetDim();
    _nNodes = _mesh->GetNumNodes();

    _nq          = quad->GetNq();
    _nSys        = _nq;
    _localNSys   = _nq / _numProcs;
    _startSysIdx = _rank * _localNSys;
    _endSysIdx   = _rank * _localNSys + _localNSys;
    if( _rank == _numProcs - 1 ) {
        _localNSys += _nq % _numProcs;
        _endSysIdx = _nSys;
    }

    _spatialOrder  = _settings->GetSpatialOrder();
    _temporalOrder = _settings->GetTemporalOrder();

    _areas                  = std::vector<double>( _nCells );
    _normals                = std::vector<double>( _nCells * _nNbr * _nDim );
    _neighbors              = std::vector<unsigned>( _nCells * _nNbr );
    _cellMidPoints          = std::vector<double>( _nCells * _nDim );
    _interfaceMidPoints     = std::vector<double>( _nCells * _nNbr * _nDim );
    _cellBoundaryTypes      = std::vector<BOUNDARY_TYPE>( _nCells );
    _relativeInterfaceMidPt = std::vector<double>( _nCells * _nNbr * _nDim );
    _relativeCellVertices   = std::vector<double>( _nCells * _nNbr * _nDim );

    // Slope
    _solDx   = std::vector<double>( _nCells * _localNSys * _nDim );
    _limiter = std::vector<double>( _nCells * _localNSys );

    // Physics
    _sigmaS = std::vector<double>( _nCells );
    _sigmaT = std::vector<double>( _nCells );
    _source = std::vector<double>( _nCells * _localNSys );

    // Quadrature
    _quadPts     = std::vector<double>( _localNSys * _nDim );
    _quadWeights = std::vector<double>( _localNSys );

    // Solution
    _sol  = std::vector<double>( _nCells * _localNSys );
    _flux = std::vector<double>( _nCells * _localNSys );

    _scalarFlux              = std::vector<double>( _nCells );
    _localMaxOrdinateOutflow = std::vector<double>( _nCells );

    auto areas           = _mesh->GetCellAreas();
    auto neighbors       = _mesh->GetNeighbours();
    auto normals         = _mesh->GetNormals();
    auto cellMidPts      = _mesh->GetCellMidPoints();
    auto interfaceMidPts = _mesh->GetInterfaceMidPoints();
    auto boundaryTypes   = _mesh->GetBoundaryTypes();
    auto nodes           = _mesh->GetNodes();
    auto cells           = _mesh->GetCells();

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        _areas[idx_cell] = areas[idx_cell];
    }

    _dT    = ComputeTimeStep( _settings->GetCFL() );
    _nIter = unsigned( _settings->GetTEnd() / _dT ) + 1;

    auto quadPoints  = quad->GetPoints();
    auto quadWeights = quad->GetWeights();

    _problem    = ProblemBase::Create( _settings, _mesh, quad );
    auto sigmaT = _problem->GetTotalXS( Vector( _nIter, 0.0 ) );
    auto sigmaS = _problem->GetScatteringXS( Vector( _nIter, 0.0 ) );
    auto source = _problem->GetExternalSource( Vector( _nIter, 0.0 ) );

    // Copy to everything to solver
    _mass = 0;
#pragma omp parallel for
    for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
        for( unsigned idx_dim = 0; idx_dim < _nDim; idx_dim++ ) {
            _quadPts[Idx2D( idx_sys, idx_dim, _nDim )] = quadPoints[idx_sys + _startSysIdx][idx_dim];
        }
        _quadWeights[idx_sys] =
            2.0 * quadWeights[idx_sys + _startSysIdx];    // Rescaling of quadweights TODO: Check if this needs general refactoring
    }

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        _cellBoundaryTypes[idx_cell] = boundaryTypes[idx_cell];

        for( unsigned idx_dim = 0; idx_dim < _nDim; idx_dim++ ) {
            _cellMidPoints[Idx2D( idx_cell, idx_dim, _nDim )] = cellMidPts[idx_cell][idx_dim];
        }

        for( unsigned idx_nbr = 0; idx_nbr < _nNbr; idx_nbr++ ) {

            _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )] = neighbors[idx_cell][idx_nbr];

            for( unsigned idx_dim = 0; idx_dim < _nDim; idx_dim++ ) {
                _normals[Idx3D( idx_cell, idx_nbr, idx_dim, _nNbr, _nDim )]            = normals[idx_cell][idx_nbr][idx_dim];
                _interfaceMidPoints[Idx3D( idx_cell, idx_nbr, idx_dim, _nNbr, _nDim )] = interfaceMidPts[idx_cell][idx_nbr][idx_dim];
                _relativeInterfaceMidPt[Idx3D( idx_cell, idx_nbr, idx_dim, _nNbr, _nDim )] =
                    _interfaceMidPoints[Idx3D( idx_cell, idx_nbr, idx_dim, _nNbr, _nDim )] - _cellMidPoints[Idx2D( idx_cell, idx_dim, _nDim )];
                _relativeCellVertices[Idx3D( idx_cell, idx_nbr, idx_dim, _nNbr, _nDim )] =
                    nodes[cells[idx_cell][idx_nbr]][idx_dim] - cellMidPts[idx_cell][idx_dim];
            }
        }

        _sigmaS[idx_cell]     = sigmaS[0][idx_cell];
        _sigmaT[idx_cell]     = sigmaT[0][idx_cell];
        _scalarFlux[idx_cell] = 0;

        for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
            _source[Idx2D( idx_cell, idx_sys, _localNSys )] = source[0][idx_cell][0];    // CAREFUL HERE hardcoded to isotropic source
            _sol[Idx2D( idx_cell, idx_sys, _localNSys )]    = 0.0;                       // initial condition is zero

            _scalarFlux[idx_cell] += _sol[Idx2D( idx_cell, idx_sys, _localNSys )] * _quadWeights[idx_sys];
        }
        // _mass += _scalarFlux[idx_cell] * _areas[idx_cell];
    }
#ifdef BUILD_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif
    SetGhostCells();
    if( _rank == 0 ) {
        PrepareScreenOutput();     // Screen Output
        PrepareHistoryOutput();    // History Output
    }
#ifdef BUILD_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif
    delete quad;

    // Initialiye QOIS
    _mass    = 0;
    _rmsFlux = 0;

    // Lattice
    {
        _curAbsorptionLattice    = 0;
        _totalAbsorptionLattice  = 0;
        _curMaxAbsorptionLattice = 0;
        _curScalarOutflow        = 0;
        _totalScalarOutflow      = 0;
        _curMaxOrdinateOutflow   = 0;
        _curScalarOutflowPeri1   = 0;
        _totalScalarOutflowPeri1 = 0;
        _curScalarOutflowPeri2   = 0;
        _totalScalarOutflowPeri2 = 0;
        ComputeCellsPerimeterLattice();
    }
    // Hohlraum
    {
        _totalAbsorptionHohlraumCenter     = 0;
        _totalAbsorptionHohlraumVertical   = 0;
        _totalAbsorptionHohlraumHorizontal = 0;
        _curAbsorptionHohlraumCenter       = 0;
        _curAbsorptionHohlraumVertical     = 0;
        _curAbsorptionHohlraumHorizontal   = 0;
        _varAbsorptionHohlraumGreen        = 0;

        unsigned n_probes = 4;
        if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) n_probes = 4;
        if( _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) n_probes = 2;
        //_probingMomentsTimeIntervals = 10;
        _probingMoments = std::vector<double>( n_probes * 3, 0. );    // 10 time horizons

        if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) {
            _probingCellsHohlraum = {
                _mesh->GetCellsofBall( -0.4, 0., 0.01 ),
                _mesh->GetCellsofBall( 0.4, 0., 0.01 ),
                _mesh->GetCellsofBall( 0., -0.5, 0.01 ),
                _mesh->GetCellsofBall( 0., 0.5, 0.01 ),
            };
        }
        else if( _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) {
            _probingCellsHohlraum = {
                _mesh->GetCellsofBall( 0.4, 0., 0.01 ),
                _mesh->GetCellsofBall( 0., 0.5, 0.01 ),
            };
        }

        // Green
        _thicknessGreen = 0.05;

        if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) {
            _centerGreen           = { _settings->GetPosXCenterGreenHohlraum(), _settings->GetPosYCenterGreenHohlraum() };
            _cornerUpperLeftGreen  = { -0.2 + _centerGreen[0], 0.4 + _centerGreen[1] };
            _cornerLowerLeftGreen  = { -0.2 + _centerGreen[0], -0.4 + _centerGreen[1] };
            _cornerUpperRightGreen = { 0.2 + _centerGreen[0], 0.4 + _centerGreen[1] };
            _cornerLowerRightGreen = { 0.2 + _centerGreen[0], -0.4 + _centerGreen[1] };
        }
        else {
            _centerGreen           = { 0.0, 0.0 };
            _cornerUpperLeftGreen  = { 0., 0.4 - _thicknessGreen / 2.0 };
            _cornerLowerLeftGreen  = { 0., +_thicknessGreen / 2.0 };
            _cornerUpperRightGreen = { 0.2 - _thicknessGreen / 2.0, 0.4 - _thicknessGreen / 2.0 };
            _cornerLowerRightGreen = { 0.2 - _thicknessGreen / 2.0, 0. + _thicknessGreen / 2.0 };
        }

        _nProbingCellsLineGreen = _settings->GetNumProbingCellsLineHohlraum();

        SetProbingCellsLineGreen();    // ONLY FOR HOHLRAUM

        _absorptionValsIntegrated    = std::vector<double>( _nProbingCellsLineGreen, 0.0 );
        _varAbsorptionValsIntegrated = std::vector<double>( _nProbingCellsLineGreen, 0.0 );
    }
}

SNSolverHPC::~SNSolverHPC() {
    delete _mesh;
    delete _problem;
}

void SNSolverHPC::Solve() {

    // --- Preprocessing ---
    if( _rank == 0 ) {
        PrepareVolumeOutput();
        DrawPreSolverOutput();
    }

    auto start = std::chrono::high_resolution_clock::now();    // Start timing

    std::chrono::duration<double> duration;
    // Loop over energies (pseudo-time of continuous slowing down approach)
    for( unsigned iter = 0; iter < _nIter; iter++ ) {
        if( iter == _nIter - 1 ) {    // last iteration
            _dT = _settings->GetTEnd() - iter * _dT;
        }
        if( _temporalOrder == 2 ) {
            std::vector<double> solRK0( _sol );
            ( _spatialOrder == 2 ) ? FluxOrder2() : FluxOrder1();
            FVMUpdate();
            ( _spatialOrder == 2 ) ? FluxOrder2() : FluxOrder1();
            FVMUpdate();
#pragma omp parallel for
            for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
#pragma omp simd
                for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
                    _sol[Idx2D( idx_cell, idx_sys, _localNSys )] =
                        0.5 * ( solRK0[Idx2D( idx_cell, idx_sys, _localNSys )] +
                                _sol[Idx2D( idx_cell, idx_sys, _localNSys )] );    // Solution averaging with HEUN
                }
            }
        }
        else {
            ( _spatialOrder == 2 ) ? FluxOrder2() : FluxOrder1();
            FVMUpdate();
        }
        IterPostprocessing();
        // --- Wall time measurement
        duration  = std::chrono::high_resolution_clock::now() - start;
        _currTime = std::chrono::duration_cast<std::chrono::duration<double>>( duration ).count();

        // --- Write Output ---
        if( _rank == 0 ) {

            WriteScalarOutput( iter );

            // --- Update Scalar Fluxes

            // --- Print Output ---
            PrintScreenOutput( iter );
            PrintHistoryOutput( iter );
            PrintVolumeOutput( iter );
        }
#ifdef BUILD_MPI
        MPI_Barrier( MPI_COMM_WORLD );
#endif
    }
    // --- Postprocessing ---
    if( _rank == 0 ) {
        DrawPostSolverOutput();
    }
}

void SNSolverHPC::PrintVolumeOutput( int idx_iter ) {
    if( _settings->GetVolumeOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) {
        WriteVolumeOutput( idx_iter );
        ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( idx_iter ), _outputFields, _outputFieldNames, _mesh );    // slow
    }
    if( idx_iter == (int)_nIter - 1 ) {    // Last iteration write without suffix.
        WriteVolumeOutput( idx_iter );
        ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh );
    }
}

void SNSolverHPC::FluxOrder2() {

    double const eps = 1e-10;

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {     // Compute Slopes
        if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NONE ) {    // skip ghost cells
#pragma omp simd
            for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {

                _limiter[Idx2D( idx_cell, idx_sys, _localNSys )] = 1.;    // limiter should be zero at boundary

                _solDx[Idx3D( idx_cell, idx_sys, 0, _localNSys, _nDim )] = 0.;
                _solDx[Idx3D( idx_cell, idx_sys, 1, _localNSys, _nDim )] = 0.;

                double solInterfaceAvg = 0.0;
                for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {    // Compute slopes and mininum and maximum
                    unsigned idx_nbr_glob = _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )];

                    // Slopes
                    solInterfaceAvg = 0.5 * ( _sol[Idx2D( idx_cell, idx_sys, _localNSys )] + _sol[Idx2D( idx_nbr_glob, idx_sys, _localNSys )] );
                    _solDx[Idx3D( idx_cell, idx_sys, 0, _localNSys, _nDim )] +=
                        solInterfaceAvg * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )];
                    _solDx[Idx3D( idx_cell, idx_sys, 1, _localNSys, _nDim )] +=
                        solInterfaceAvg * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];
                }

                _solDx[Idx3D( idx_cell, idx_sys, 0, _localNSys, _nDim )] /= _areas[idx_cell];
                _solDx[Idx3D( idx_cell, idx_sys, 1, _localNSys, _nDim )] /= _areas[idx_cell];
            }
        }
    }
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {     // Compute Limiter
        if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NONE ) {    // skip ghost cells

#pragma omp simd
            for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {

                double gaussPoint = 0;
                double r          = 0;
                double minSol     = _sol[Idx2D( idx_cell, idx_sys, _localNSys )];
                double maxSol     = _sol[Idx2D( idx_cell, idx_sys, _localNSys )];

                for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {    // Compute slopes and mininum and maximum
                    unsigned idx_nbr_glob = _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )];
                    // Compute ptswise max and minimum solultion values of current and neighbor cells
                    maxSol = std::max( _sol[Idx2D( idx_nbr_glob, idx_sys, _localNSys )], maxSol );
                    minSol = std::min( _sol[Idx2D( idx_nbr_glob, idx_sys, _localNSys )], minSol );
                }

                for( unsigned idx_nbr = 0; idx_nbr < _nNbr; idx_nbr++ ) {    // Compute limiter, see https://arxiv.org/pdf/1710.07187.pdf

                    // Compute test value at cell vertex, called gaussPt
                    gaussPoint =
                        _solDx[Idx3D( idx_cell, idx_sys, 0, _localNSys, _nDim )] *
                            _relativeCellVertices[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                        _solDx[Idx3D( idx_cell, idx_sys, 1, _localNSys, _nDim )] * _relativeCellVertices[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];

                    // BARTH-JESPERSEN LIMITER
                    // r = ( gaussPoint > 0 ) ? std::min( ( maxSol - _sol[Idx2D( idx_cell, idx_sys, _localNSys )] ) / ( gaussPoint + eps ), 1.0 )
                    //                       : std::min( ( minSol - _sol[Idx2D( idx_cell, idx_sys, _localNSys )] ) / ( gaussPoint - eps ), 1.0 );
                    //
                    // r = ( std::abs( gaussPoint ) < eps ) ? 1 : r;
                    //_limiter[Idx2D( idx_cell, idx_sys, _localNSys )] = std::min( r, _limiter[Idx2D( idx_cell, idx_sys, _localNSys )] );

                    // VENKATAKRISHNAN LIMITER
                    double delta1Max = maxSol - _sol[Idx2D( idx_cell, idx_sys, _localNSys )];
                    double delta1Min = minSol - _sol[Idx2D( idx_cell, idx_sys, _localNSys )];

                    r = ( gaussPoint > 0 ) ? ( ( delta1Max * delta1Max + _areas[idx_cell] ) * gaussPoint + 2 * gaussPoint * gaussPoint * delta1Max ) /
                                                 ( delta1Max * delta1Max + 2 * gaussPoint * gaussPoint + delta1Max * gaussPoint + _areas[idx_cell] ) /
                                                 ( gaussPoint + eps )
                                           : ( ( delta1Min * delta1Min + _areas[idx_cell] ) * gaussPoint + 2 * gaussPoint * gaussPoint * delta1Min ) /
                                                 ( delta1Min * delta1Min + 2 * gaussPoint * gaussPoint + delta1Min * gaussPoint + _areas[idx_cell] ) /
                                                 ( gaussPoint - eps );

                    r = ( std::abs( gaussPoint ) < eps ) ? 1 : r;

                    _limiter[Idx2D( idx_cell, idx_sys, _localNSys )] = std::min( r, _limiter[Idx2D( idx_cell, idx_sys, _localNSys )] );
                }
            }
        }
        else {
#pragma omp simd
            for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
                _limiter[Idx2D( idx_cell, idx_sys, _localNSys )]         = 0.;    // limiter should be zero at boundary
                _solDx[Idx3D( idx_cell, idx_sys, 0, _localNSys, _nDim )] = 0.;
                _solDx[Idx3D( idx_cell, idx_sys, 1, _localNSys, _nDim )] = 0.;
            }
        }
    }

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {    // Compute Flux
#pragma omp simd
        for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
            _flux[Idx2D( idx_cell, idx_sys, _localNSys )] = 0.;
        }

        // Fluxes
        for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {
            if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )] == _nCells ) {
#pragma omp simd
                for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
                    double localInner = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                        _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];
                    if( localInner > 0 ) {
                        _flux[Idx2D( idx_cell, idx_sys, _localNSys )] += localInner * _sol[Idx2D( idx_cell, idx_sys, _localNSys )];
                    }
                    else {
                        double ghostCellValue = _ghostCells[idx_cell][idx_sys];    // fixed boundary
                        _flux[Idx2D( idx_cell, idx_sys, _localNSys )] += localInner * ghostCellValue;
                    }
                }
            }
            else {
// Second order
#pragma omp simd
                for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
                    // store flux contribution on psiNew_sigmaS to save memory
                    unsigned idx_nbr_glob = _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )];    // global idx of neighbor cell
                    double localInner     = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                        _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];

                    _flux[Idx2D( idx_cell, idx_sys, _localNSys )] +=
                        ( localInner > 0 ) ? localInner * ( _sol[Idx2D( idx_cell, idx_sys, _localNSys )] +
                                                            _limiter[Idx2D( idx_cell, idx_sys, _localNSys )] *
                                                                ( _solDx[Idx3D( idx_cell, idx_sys, 0, _localNSys, _nDim )] *
                                                                      _relativeInterfaceMidPt[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                                                  _solDx[Idx3D( idx_cell, idx_sys, 1, _localNSys, _nDim )] *
                                                                      _relativeInterfaceMidPt[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )] ) )
                                           : localInner * ( _sol[Idx2D( idx_nbr_glob, idx_sys, _localNSys )] +
                                                            _limiter[Idx2D( idx_nbr_glob, idx_sys, _localNSys )] *
                                                                ( _solDx[Idx3D( idx_nbr_glob, idx_sys, 0, _localNSys, _nDim )] *
                                                                      ( _interfaceMidPoints[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] -
                                                                        _cellMidPoints[Idx2D( idx_nbr_glob, 0, _nDim )] ) +
                                                                  _solDx[Idx3D( idx_nbr_glob, idx_sys, 1, _localNSys, _nDim )] *
                                                                      ( _interfaceMidPoints[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )] -
                                                                        _cellMidPoints[Idx2D( idx_nbr_glob, 1, _nDim )] ) ) );
                }
            }
        }
    }
}

void SNSolverHPC::FluxOrder1() {

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {

#pragma omp simd
        for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
            _flux[Idx2D( idx_cell, idx_sys, _localNSys )] = 0.0;    // Reset temporary variable
        }

        // Fluxes
        for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {
            if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )] == _nCells ) {
#pragma omp simd
                for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
                    double localInner = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                        _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];
                    if( localInner > 0 ) {
                        _flux[Idx2D( idx_cell, idx_sys, _localNSys )] += localInner * _sol[Idx2D( idx_cell, idx_sys, _localNSys )];
                    }
                    else {
                        double ghostCellValue = _ghostCells[idx_cell][idx_sys];    // fixed boundary
                        _flux[Idx2D( idx_cell, idx_sys, _localNSys )] += localInner * ghostCellValue;
                    }
                }
            }
            else {
                unsigned nbr_glob = _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )];    // global idx of neighbor cell
#pragma omp simd
                for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {

                    double localInner = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                        _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];

                    _flux[Idx2D( idx_cell, idx_sys, _localNSys )] += ( localInner > 0 ) ? localInner * _sol[Idx2D( idx_cell, idx_sys, _localNSys )]
                                                                                        : localInner * _sol[Idx2D( nbr_glob, idx_sys, _localNSys )];
                }
            }
        }
    }
}

void SNSolverHPC::FVMUpdate() {
    _mass    = 0.0;
    _rmsFlux = 0.0;
    std::vector<double> temp_scalarFlux( _nCells );    // for MPI allreduce

#pragma omp parallel for reduction( + : _mass, _rmsFlux )
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {

#pragma omp simd
        for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
            // Update
            _sol[Idx2D( idx_cell, idx_sys, _localNSys )] =
                ( 1 - _dT * _sigmaT[idx_cell] ) * _sol[Idx2D( idx_cell, idx_sys, _localNSys )] -
                _dT / _areas[idx_cell] * _flux[Idx2D( idx_cell, idx_sys, _localNSys )] +
                _dT * ( _sigmaS[idx_cell] * _scalarFlux[idx_cell] / ( 2 * M_PI ) + _source[Idx2D( idx_cell, idx_sys, _localNSys )] );
        }
        double localScalarFlux = 0;

#pragma omp simd reduction( + : localScalarFlux )
        for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
            _sol[Idx2D( idx_cell, idx_sys, _localNSys )] = std::max( _sol[Idx2D( idx_cell, idx_sys, _localNSys )], 0.0 );
            localScalarFlux += _sol[Idx2D( idx_cell, idx_sys, _localNSys )] * _quadWeights[idx_sys];
        }
        _mass += localScalarFlux * _areas[idx_cell];
        _rmsFlux += ( localScalarFlux - _scalarFlux[idx_cell] ) * ( localScalarFlux - _scalarFlux[idx_cell] );
        temp_scalarFlux[idx_cell] = localScalarFlux;    // set flux
    }
// MPI Allreduce: _scalarFlux
#ifdef BUILD_MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Allreduce( temp_scalarFlux.data(), _scalarFlux.data(), _nCells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );
#endif
#ifndef BUILD_MPI
    _scalarFlux = temp_scalarFlux;
#endif
}

void SNSolverHPC::IterPostprocessing() {
    // ALREDUCE NEEDED

    _curAbsorptionLattice            = 0.0;
    _curScalarOutflow                = 0.0;
    _curScalarOutflowPeri1           = 0.0;
    _curScalarOutflowPeri2           = 0.0;
    _curAbsorptionHohlraumCenter     = 0.0;    // Green and blue areas of symmetric hohlraum
    _curAbsorptionHohlraumVertical   = 0.0;    // Red areas of symmetric hohlraum
    _curAbsorptionHohlraumHorizontal = 0.0;    // Black areas of symmetric hohlraum
    _varAbsorptionHohlraumGreen      = 0.0;
    double a_g                       = 0.0;

#pragma omp parallel for reduction( + : _curAbsorptionLattice,                                                                                       \
                                        _curScalarOutflow,                                                                                           \
                                        _curScalarOutflowPeri1,                                                                                      \
                                        _curScalarOutflowPeri2,                                                                                      \
                                        _curAbsorptionHohlraumCenter,                                                                                \
                                        _curAbsorptionHohlraumVertical,                                                                              \
                                        _curAbsorptionHohlraumHorizontal,                                                                            \
                                        a_g ) reduction( max : _curMaxOrdinateOutflow, _curMaxAbsorptionLattice )
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {

        if( _settings->GetProblemName() == PROBLEM_Lattice || _settings->GetProblemName() == PROBLEM_HalfLattice ) {
            if( IsAbsorptionLattice( _cellMidPoints[Idx2D( idx_cell, 0, _nDim )], _cellMidPoints[Idx2D( idx_cell, 1, _nDim )] ) ) {
                double sigmaAPsi = _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * _areas[idx_cell];
                _curAbsorptionLattice += sigmaAPsi;
                _curMaxAbsorptionLattice = ( _curMaxAbsorptionLattice < sigmaAPsi ) ? sigmaAPsi : _curMaxAbsorptionLattice;
            }
        }

        if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum || _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) {

            double x = _cellMidPoints[Idx2D( idx_cell, 0, _nDim )];
            double y = _cellMidPoints[Idx2D( idx_cell, 1, _nDim )];
            _curAbsorptionLattice += _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * _areas[idx_cell];
            if( x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                y > -0.4 + _settings->GetPosYCenterGreenHohlraum() && y < 0.4 + _settings->GetPosYCenterGreenHohlraum() ) {
                _curAbsorptionHohlraumCenter += _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * _areas[idx_cell];
            }
            if( ( x < _settings->GetPosRedLeftBorderHohlraum() && y > _settings->GetPosRedLeftBottomHohlraum() &&
                  y < _settings->GetPosRedLeftTopHohlraum() ) ||
                ( x > _settings->GetPosRedRightBorderHohlraum() && y > _settings->GetPosRedLeftBottomHohlraum() &&
                  y < _settings->GetPosRedRightTopHohlraum() ) ) {
                _curAbsorptionHohlraumVertical += _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * _areas[idx_cell];
            }
            if( y > 0.6 || y < -0.6 ) {
                _curAbsorptionHohlraumHorizontal += _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * _areas[idx_cell];
            }

            // Variation in absorption of center (part 1)
            // green area 1 (lower boundary)
            bool green1 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                          y > -0.4 + _settings->GetPosYCenterGreenHohlraum() && y < -0.35 + _settings->GetPosYCenterGreenHohlraum();
            // green area 2 (upper boundary)
            bool green2 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                          y > 0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.4 + _settings->GetPosYCenterGreenHohlraum();
            // green area 3 (left boundary)
            bool green3 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < -0.15 + _settings->GetPosXCenterGreenHohlraum() &&
                          y > -0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.35 + _settings->GetPosYCenterGreenHohlraum();
            // green area 4 (right boundary)
            bool green4 = x > 0.15 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                          y > -0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.35 + _settings->GetPosYCenterGreenHohlraum();
            if( green1 || green2 || green3 || green4 ) {
                a_g += ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * _scalarFlux[idx_cell] * _areas[idx_cell];
            }
        }

        if( _settings->GetProblemName() == PROBLEM_Lattice ) {
            // Outflow out of inner and middle perimeter
            if( _isPerimeterLatticeCell1[idx_cell] ) {    // inner
                for( unsigned idx_nbr_helper = 0; idx_nbr_helper < _cellsLatticePerimeter1[idx_cell].size(); ++idx_nbr_helper ) {
#pragma omp simd reduction( + : _curScalarOutflowPeri1 )
                    for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
                        double localInner = _quadPts[Idx2D( idx_sys, 0, _nDim )] *
                                                _normals[Idx3D( idx_cell, _cellsLatticePerimeter1[idx_cell][idx_nbr_helper], 0, _nNbr, _nDim )] +
                                            _quadPts[Idx2D( idx_sys, 1, _nDim )] *
                                                _normals[Idx3D( idx_cell, _cellsLatticePerimeter1[idx_cell][idx_nbr_helper], 1, _nNbr, _nDim )];
                        // Find outward facing transport directions

                        if( localInner > 0.0 ) {
                            _curScalarOutflowPeri1 +=
                                localInner * _sol[Idx2D( idx_cell, idx_sys, _localNSys )] * _quadWeights[idx_sys];    // Integrate flux
                        }
                    }
                }
            }
            if( _isPerimeterLatticeCell2[idx_cell] ) {    // middle
                for( unsigned idx_nbr_helper = 0; idx_nbr_helper < _cellsLatticePerimeter2[idx_cell].size(); ++idx_nbr_helper ) {
#pragma omp simd reduction( + : _curScalarOutflowPeri2 )
                    for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
                        double localInner = _quadPts[Idx2D( idx_sys, 0, _nDim )] *
                                                _normals[Idx3D( idx_cell, _cellsLatticePerimeter2[idx_cell][idx_nbr_helper], 0, _nNbr, _nDim )] +
                                            _quadPts[Idx2D( idx_sys, 1, _nDim )] *
                                                _normals[Idx3D( idx_cell, _cellsLatticePerimeter2[idx_cell][idx_nbr_helper], 1, _nNbr, _nDim )];
                        // Find outward facing transport directions

                        if( localInner > 0.0 ) {
                            _curScalarOutflowPeri2 +=
                                localInner * _sol[Idx2D( idx_cell, idx_sys, _localNSys )] * _quadWeights[idx_sys];    // Integrate flux
                        }
                    }
                }
            }
        }
        // Outflow out of domain
        if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NEUMANN ) {
            // Iterate over face cell faces
            double currOrdinatewiseOutflow = 0.0;

            for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {
                // Find face that points outward
                if( _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )] == _nCells ) {
#pragma omp simd reduction( + : _curScalarOutflow ) reduction( max : _curMaxOrdinateOutflow )
                    for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {

                        double localInner = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                            _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];
                        // Find outward facing transport directions

                        if( localInner > 0.0 ) {
                            _curScalarOutflow +=
                                localInner * _sol[Idx2D( idx_cell, idx_sys, _localNSys )] * _quadWeights[idx_sys];    // Integrate flux

                            currOrdinatewiseOutflow =
                                _sol[Idx2D( idx_cell, idx_sys, _localNSys )] * localInner /
                                sqrt( (
                                    _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                    _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )] ) );

                            _curMaxOrdinateOutflow =
                                ( currOrdinatewiseOutflow > _curMaxOrdinateOutflow ) ? currOrdinatewiseOutflow : _curMaxOrdinateOutflow;
                        }
                    }
                }
            }
        }
    }
// MPI Allreduce
#ifdef BUILD_MPI
    double tmp_curScalarOutflow      = 0.0;
    double tmp_curScalarOutflowPeri1 = 0.0;
    double tmp_curScalarOutflowPeri2 = 0.0;
    double tmp_mass                  = 0.0;
    double tmp_rmsFlux               = 0.0;
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Allreduce( &_curScalarOutflow, &tmp_curScalarOutflow, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    _curScalarOutflow = tmp_curScalarOutflow;
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Allreduce( &_curScalarOutflowPeri1, &tmp_curScalarOutflowPeri1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    _curScalarOutflowPeri1 = tmp_curScalarOutflowPeri1;
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Allreduce( &_curScalarOutflowPeri2, &tmp_curScalarOutflowPeri2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    _curScalarOutflowPeri2 = tmp_curScalarOutflowPeri2;
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Allreduce( &_mass, &tmp_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    _mass = tmp_mass;
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Allreduce( &_rmsFlux, &tmp_rmsFlux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    _rmsFlux = tmp_rmsFlux;
    MPI_Barrier( MPI_COMM_WORLD );
#endif
    // Variation absorption (part II)
    if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) {
        unsigned n_probes = 4;
        std::vector<double> temp_probingMoments( 3 * n_probes );    // for MPI allreduce

#pragma omp parallel for reduction( + : _varAbsorptionHohlraumGreen )
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
            double x = _cellMidPoints[Idx2D( idx_cell, 0, _nDim )];
            double y = _cellMidPoints[Idx2D( idx_cell, 1, _nDim )];

            // green area 1 (lower boundary)
            bool green1 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                          y > -0.4 + _settings->GetPosYCenterGreenHohlraum() && y < -0.35 + _settings->GetPosYCenterGreenHohlraum();
            // green area 2 (upper boundary)
            bool green2 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                          y > 0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.4 + _settings->GetPosYCenterGreenHohlraum();
            // green area 3 (left boundary)
            bool green3 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < -0.15 + _settings->GetPosXCenterGreenHohlraum() &&
                          y > -0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.35 + _settings->GetPosYCenterGreenHohlraum();
            // green area 4 (right boundary)
            bool green4 = x > 0.15 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                          y > -0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.35 + _settings->GetPosYCenterGreenHohlraum();
            if( green1 || green2 || green3 || green4 ) {
                _varAbsorptionHohlraumGreen += ( a_g - _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) ) *
                                               ( a_g - _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) ) * _areas[idx_cell];
            }
        }
        // Probes value moments
        // #pragma omp parallel for
        for( unsigned idx_probe = 0; idx_probe < n_probes; idx_probe++ ) {    // Loop over probing cells
            temp_probingMoments[Idx2D( idx_probe, 0, 3 )] = 0.0;
            temp_probingMoments[Idx2D( idx_probe, 1, 3 )] = 0.0;
            temp_probingMoments[Idx2D( idx_probe, 2, 3 )] = 0.0;
            // for( unsigned idx_ball = 0; idx_ball < _probingCellsHohlraum[idx_probe].size(); idx_ball++ ) {
            //   std::cout << idx_ball << _areas[_probingCellsHohlraum[idx_probe][idx_ball]] / ( 0.01 * 0.01 * M_PI ) << std::endl;
            //}
            for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
                for( unsigned idx_ball = 0; idx_ball < _probingCellsHohlraum[idx_probe].size(); idx_ball++ ) {
                    temp_probingMoments[Idx2D( idx_probe, 0, 3 )] += _sol[Idx2D( _probingCellsHohlraum[idx_probe][idx_ball], idx_sys, _localNSys )] *
                                                                     _quadWeights[idx_sys] * _areas[_probingCellsHohlraum[idx_probe][idx_ball]] /
                                                                     ( 0.01 * 0.01 * M_PI );
                    temp_probingMoments[Idx2D( idx_probe, 1, 3 )] +=
                        _quadPts[Idx2D( idx_sys, 0, _nDim )] * _sol[Idx2D( _probingCellsHohlraum[idx_probe][idx_ball], idx_sys, _localNSys )] *
                        _quadWeights[idx_sys] * _areas[_probingCellsHohlraum[idx_probe][idx_ball]] / ( 0.01 * 0.01 * M_PI );
                    temp_probingMoments[Idx2D( idx_probe, 2, 3 )] +=
                        _quadPts[Idx2D( idx_sys, 1, _nDim )] * _sol[Idx2D( _probingCellsHohlraum[idx_probe][idx_ball], idx_sys, _localNSys )] *
                        _quadWeights[idx_sys] * _areas[_probingCellsHohlraum[idx_probe][idx_ball]] / ( 0.01 * 0.01 * M_PI );
                }
            }
        }

        // probe values green
        ComputeQOIsGreenProbingLine();
#ifdef BUILD_MPI
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Allreduce( temp_probingMoments.data(), _probingMoments.data(), 3 * n_probes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        MPI_Barrier( MPI_COMM_WORLD );
#endif
#ifndef BUILD_MPI
        for( unsigned idx_probe = 0; idx_probe < n_probes; idx_probe++ ) {    // Loop over probing cells
            _probingMoments[Idx2D( idx_probe, 0, 3 )] = temp_probingMoments[Idx2D( idx_probe, 0, 3 )];
            _probingMoments[Idx2D( idx_probe, 1, 3 )] = temp_probingMoments[Idx2D( idx_probe, 1, 3 )];
            _probingMoments[Idx2D( idx_probe, 2, 3 )] = temp_probingMoments[Idx2D( idx_probe, 2, 3 )];
        }
#endif
    }

    // Update time integral values on rank 0
    if( _rank == 0 ) {
        _totalScalarOutflow += _curScalarOutflow * _dT;
        _totalScalarOutflowPeri1 += _curScalarOutflowPeri1 * _dT;
        _totalScalarOutflowPeri2 += _curScalarOutflowPeri2 * _dT;
        _totalAbsorptionLattice += _curAbsorptionLattice * _dT;

        _totalAbsorptionHohlraumCenter += _curAbsorptionHohlraumCenter * _dT;
        _totalAbsorptionHohlraumVertical += _curAbsorptionHohlraumVertical * _dT;
        _totalAbsorptionHohlraumHorizontal += _curAbsorptionHohlraumHorizontal * _dT;

        _rmsFlux = sqrt( _rmsFlux );
    }
}

bool SNSolverHPC::IsAbsorptionLattice( double x, double y ) const {
    // Check whether pos is inside absorbing squares
    double xy_corrector = -3.5;
    std::vector<double> lbounds{ 1 + xy_corrector, 2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector };
    std::vector<double> ubounds{ 2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector, 6 + xy_corrector };
    for( unsigned k = 0; k < lbounds.size(); ++k ) {
        for( unsigned l = 0; l < lbounds.size(); ++l ) {
            if( ( l + k ) % 2 == 1 || ( k == 2 && l == 2 ) || ( k == 2 && l == 4 ) ) continue;
            if( x >= lbounds[k] && x <= ubounds[k] && y >= lbounds[l] && y <= ubounds[l] ) {
                return true;
            }
        }
    }
    return false;
}

// --- Helper ---
double SNSolverHPC::ComputeTimeStep( double cfl ) const {
    // for pseudo 1D, set timestep to dx
    double dx, dy;
    switch( _settings->GetProblemName() ) {
        case PROBLEM_Checkerboard1D:
            dx = 7.0 / (double)_nCells;
            dy = 0.3;
            return cfl * ( dx * dy ) / ( dx + dy );
            break;
        case PROBLEM_Linesource1D:     // Fallthrough
        case PROBLEM_Meltingcube1D:    // Fallthrough
        case PROBLEM_Aircavity1D:
            dx = 3.0 / (double)_nCells;
            dy = 0.3;
            return cfl * ( dx * dy ) / ( dx + dy );
            break;
        default: break;    // 2d as normal
    }
    // 2D case
    double charSize = __DBL_MAX__;    // minimum char size of all mesh cells in the mesh

#pragma omp parallel for reduction( min : charSize )
    for( unsigned j = 0; j < _nCells; j++ ) {
        double currCharSize = sqrt( _areas[j] );
        if( currCharSize < charSize ) {
            charSize = currCharSize;
        }
    }
    if( _rank == 0 ) {
        auto log         = spdlog::get( "event" );
        std::string line = "| Smallest characteristic length of a grid cell in this mesh: " + std::to_string( charSize );
        log->info( line );
        line = "| Corresponding maximal time-step: " + std::to_string( cfl * charSize );
        log->info( line );
    }
    return cfl * charSize;
}

// --- IO ----
void SNSolverHPC::PrepareScreenOutput() {
    unsigned nFields = (unsigned)_settings->GetNScreenOutput();

    _screenOutputFieldNames.resize( nFields );
    _screenOutputFields.resize( nFields );

    // Prepare all output Fields ==> Specified in option SCREEN_OUTPUT
    for( unsigned idx_field = 0; idx_field < nFields; idx_field++ ) {
        // Prepare all Output Fields per group

        // Different procedure, depending on the Group...
        switch( _settings->GetScreenOutput()[idx_field] ) {
            case MASS: _screenOutputFieldNames[idx_field] = "Mass"; break;
            case ITER: _screenOutputFieldNames[idx_field] = "Iter"; break;
            case WALL_TIME: _screenOutputFieldNames[idx_field] = "Wall time [s]"; break;
            case RMS_FLUX: _screenOutputFieldNames[idx_field] = "RMS flux"; break;
            case VTK_OUTPUT: _screenOutputFieldNames[idx_field] = "VTK out"; break;
            case CSV_OUTPUT: _screenOutputFieldNames[idx_field] = "CSV out"; break;
            case CUR_OUTFLOW: _screenOutputFieldNames[idx_field] = "Cur. outflow"; break;
            case TOTAL_OUTFLOW: _screenOutputFieldNames[idx_field] = "Tot. outflow"; break;
            case CUR_OUTFLOW_P1: _screenOutputFieldNames[idx_field] = "Cur. outflow P1"; break;
            case TOTAL_OUTFLOW_P1: _screenOutputFieldNames[idx_field] = "Tot. outflow P1"; break;
            case CUR_OUTFLOW_P2: _screenOutputFieldNames[idx_field] = "Cur. outflow P2"; break;
            case TOTAL_OUTFLOW_P2: _screenOutputFieldNames[idx_field] = "Tot. outflow P2"; break;
            case MAX_OUTFLOW: _screenOutputFieldNames[idx_field] = "Max outflow"; break;
            case CUR_PARTICLE_ABSORPTION: _screenOutputFieldNames[idx_field] = "Cur. absorption"; break;
            case TOTAL_PARTICLE_ABSORPTION: _screenOutputFieldNames[idx_field] = "Tot. absorption"; break;
            case MAX_PARTICLE_ABSORPTION: _screenOutputFieldNames[idx_field] = "Max absorption"; break;
            case TOTAL_PARTICLE_ABSORPTION_CENTER: _screenOutputFieldNames[idx_field] = "Tot. abs. center"; break;
            case TOTAL_PARTICLE_ABSORPTION_VERTICAL: _screenOutputFieldNames[idx_field] = "Tot. abs. vertical wall"; break;
            case TOTAL_PARTICLE_ABSORPTION_HORIZONTAL: _screenOutputFieldNames[idx_field] = "Tot. abs. horizontal wall"; break;
            case PROBE_MOMENT_TIME_TRACE:
                _screenOutputFieldNames[idx_field] = "Probe 1 u_0";
                idx_field++;
                _screenOutputFieldNames[idx_field] = "Probe 2 u_0";
                if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) {
                    idx_field++;
                    _screenOutputFieldNames[idx_field] = "Probe 3 u_0";
                    idx_field++;
                    _screenOutputFieldNames[idx_field] = "Probe 4 u_0";
                }
                break;
            case VAR_ABSORPTION_GREEN: _screenOutputFieldNames[idx_field] = "Var. absorption green"; break;
            default: ErrorMessages::Error( "Screen output field not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void SNSolverHPC::WriteScalarOutput( unsigned idx_iter ) {
    unsigned n_probes = 4;

    unsigned nFields                  = (unsigned)_settings->GetNScreenOutput();
    const VectorVector probingMoments = _problem->GetCurrentProbeMoment();
    // -- Screen Output
    for( unsigned idx_field = 0; idx_field < nFields; idx_field++ ) {
        // Prepare all Output Fields per group
        // Different procedure, depending on the Group...
        switch( _settings->GetScreenOutput()[idx_field] ) {
            case MASS: _screenOutputFields[idx_field] = _mass; break;
            case ITER: _screenOutputFields[idx_field] = idx_iter; break;
            case WALL_TIME: _screenOutputFields[idx_field] = _currTime; break;
            case RMS_FLUX: _screenOutputFields[idx_field] = _rmsFlux; break;
            case VTK_OUTPUT:
                _screenOutputFields[idx_field] = 0;
                if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
                    ( idx_iter == _nIter - 1 ) /* need sol at last iteration */ ) {
                    _screenOutputFields[idx_field] = 1;
                }
                break;
            case CSV_OUTPUT:
                _screenOutputFields[idx_field] = 0;
                if( ( _settings->GetHistoryOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetHistoryOutputFrequency() == 0 ) ||
                    ( idx_iter == _nIter - 1 ) /* need sol at last iteration */ ) {
                    _screenOutputFields[idx_field] = 1;
                }
                break;
            case CUR_OUTFLOW: _screenOutputFields[idx_field] = _curScalarOutflow; break;
            case TOTAL_OUTFLOW: _screenOutputFields[idx_field] = _totalScalarOutflow; break;
            case CUR_OUTFLOW_P1: _screenOutputFields[idx_field] = _curScalarOutflowPeri1; break;
            case TOTAL_OUTFLOW_P1: _screenOutputFields[idx_field] = _totalScalarOutflowPeri1; break;
            case CUR_OUTFLOW_P2: _screenOutputFields[idx_field] = _curScalarOutflowPeri2; break;
            case TOTAL_OUTFLOW_P2: _screenOutputFields[idx_field] = _totalScalarOutflowPeri2; break;
            case MAX_OUTFLOW: _screenOutputFields[idx_field] = _curMaxOrdinateOutflow; break;
            case CUR_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _curAbsorptionLattice; break;
            case TOTAL_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _totalAbsorptionLattice; break;
            case MAX_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _curMaxAbsorptionLattice; break;
            case TOTAL_PARTICLE_ABSORPTION_CENTER: _screenOutputFields[idx_field] = _totalAbsorptionHohlraumCenter; break;
            case TOTAL_PARTICLE_ABSORPTION_VERTICAL: _screenOutputFields[idx_field] = _totalAbsorptionHohlraumVertical; break;
            case TOTAL_PARTICLE_ABSORPTION_HORIZONTAL: _screenOutputFields[idx_field] = _totalAbsorptionHohlraumHorizontal; break;
            case PROBE_MOMENT_TIME_TRACE:
                if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) n_probes = 4;
                if( _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) n_probes = 2;
                for( unsigned i = 0; i < n_probes; i++ ) {
                    _screenOutputFields[idx_field] = _probingMoments[Idx2D( i, 0, 3 )];
                    idx_field++;
                }
                idx_field--;
                break;
            case VAR_ABSORPTION_GREEN: _screenOutputFields[idx_field] = _varAbsorptionHohlraumGreen; break;
            default: ErrorMessages::Error( "Screen output group not defined!", CURRENT_FUNCTION ); break;
        }
    }

    // --- History output ---
    nFields = (unsigned)_settings->GetNHistoryOutput();

    std::vector<SCALAR_OUTPUT> screenOutputFields = _settings->GetScreenOutput();
    for( unsigned idx_field = 0; idx_field < nFields; idx_field++ ) {

        // Prepare all Output Fields per group
        // Different procedure, depending on the Group...
        switch( _settings->GetHistoryOutput()[idx_field] ) {
            case MASS: _historyOutputFields[idx_field] = _mass; break;
            case ITER: _historyOutputFields[idx_field] = idx_iter; break;
            case WALL_TIME: _historyOutputFields[idx_field] = _currTime; break;
            case RMS_FLUX: _historyOutputFields[idx_field] = _rmsFlux; break;
            case VTK_OUTPUT:
                _historyOutputFields[idx_field] = 0;
                if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
                    ( idx_iter == _nIter - 1 ) /* need sol at last iteration */ ) {
                    _historyOutputFields[idx_field] = 1;
                }
                break;

            case CSV_OUTPUT:
                _historyOutputFields[idx_field] = 0;
                if( ( _settings->GetHistoryOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetHistoryOutputFrequency() == 0 ) ||
                    ( idx_iter == _nIter - 1 ) /* need sol at last iteration */ ) {
                    _historyOutputFields[idx_field] = 1;
                }
                break;
            case CUR_OUTFLOW: _historyOutputFields[idx_field] = _curScalarOutflow; break;
            case TOTAL_OUTFLOW: _historyOutputFields[idx_field] = _totalScalarOutflow; break;
            case CUR_OUTFLOW_P1: _historyOutputFields[idx_field] = _curScalarOutflowPeri1; break;
            case TOTAL_OUTFLOW_P1: _historyOutputFields[idx_field] = _totalScalarOutflowPeri1; break;
            case CUR_OUTFLOW_P2: _historyOutputFields[idx_field] = _curScalarOutflowPeri2; break;
            case TOTAL_OUTFLOW_P2: _historyOutputFields[idx_field] = _totalScalarOutflowPeri2; break;
            case MAX_OUTFLOW: _historyOutputFields[idx_field] = _curMaxOrdinateOutflow; break;
            case CUR_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _curAbsorptionLattice; break;
            case TOTAL_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _totalAbsorptionLattice; break;
            case MAX_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _curMaxAbsorptionLattice; break;
            case TOTAL_PARTICLE_ABSORPTION_CENTER: _historyOutputFields[idx_field] = _totalAbsorptionHohlraumCenter; break;
            case TOTAL_PARTICLE_ABSORPTION_VERTICAL: _historyOutputFields[idx_field] = _totalAbsorptionHohlraumVertical; break;
            case TOTAL_PARTICLE_ABSORPTION_HORIZONTAL: _historyOutputFields[idx_field] = _totalAbsorptionHohlraumHorizontal; break;
            case PROBE_MOMENT_TIME_TRACE:
                if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) n_probes = 4;
                if( _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) n_probes = 2;
                for( unsigned i = 0; i < n_probes; i++ ) {
                    for( unsigned j = 0; j < 3; j++ ) {
                        _historyOutputFields[idx_field] = _probingMoments[Idx2D( i, j, 3 )];
                        idx_field++;
                    }
                }
                idx_field--;
                break;
            case VAR_ABSORPTION_GREEN: _historyOutputFields[idx_field] = _varAbsorptionHohlraumGreen; break;
            case VAR_ABSORPTION_GREEN_LINE:
                for( unsigned i = 0; i < _settings->GetNumProbingCellsLineHohlraum(); i++ ) {
                    _historyOutputFieldNames[idx_field] = _varAbsorptionValsIntegrated[i];
                    idx_field++;
                }
                idx_field--;
                break;
            default: ErrorMessages::Error( "History output group not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void SNSolverHPC::PrintScreenOutput( unsigned idx_iter ) {
    auto log = spdlog::get( "event" );

    unsigned strLen  = 15;    // max width of one column
    char paddingChar = ' ';

    // assemble the line to print
    std::string lineToPrint = "| ";
    std::string tmp;
    for( unsigned idx_field = 0; idx_field < _settings->GetNScreenOutput(); idx_field++ ) {
        tmp = std::to_string( _screenOutputFields[idx_field] );

        // Format outputs correctly
        std::vector<SCALAR_OUTPUT> integerFields    = { ITER };
        std::vector<SCALAR_OUTPUT> scientificFields = { RMS_FLUX,
                                                        MASS,
                                                        WALL_TIME,
                                                        CUR_OUTFLOW,
                                                        TOTAL_OUTFLOW,
                                                        CUR_OUTFLOW_P1,
                                                        TOTAL_OUTFLOW_P1,
                                                        CUR_OUTFLOW_P2,
                                                        TOTAL_OUTFLOW_P2,
                                                        MAX_OUTFLOW,
                                                        CUR_PARTICLE_ABSORPTION,
                                                        TOTAL_PARTICLE_ABSORPTION,
                                                        MAX_PARTICLE_ABSORPTION,
                                                        TOTAL_PARTICLE_ABSORPTION_CENTER,
                                                        TOTAL_PARTICLE_ABSORPTION_VERTICAL,
                                                        TOTAL_PARTICLE_ABSORPTION_HORIZONTAL,
                                                        PROBE_MOMENT_TIME_TRACE,
                                                        VAR_ABSORPTION_GREEN,
                                                        VAR_ABSORPTION_GREEN_LINE };
        std::vector<SCALAR_OUTPUT> booleanFields    = { VTK_OUTPUT, CSV_OUTPUT };

        if( !( integerFields.end() == std::find( integerFields.begin(), integerFields.end(), _settings->GetScreenOutput()[idx_field] ) ) ) {
            tmp = std::to_string( (int)_screenOutputFields[idx_field] );
        }
        else if( !( booleanFields.end() == std::find( booleanFields.begin(), booleanFields.end(), _settings->GetScreenOutput()[idx_field] ) ) ) {
            tmp = "no";
            if( (bool)_screenOutputFields[idx_field] ) tmp = "yes";
        }
        else if( !( scientificFields.end() ==
                    std::find( scientificFields.begin(), scientificFields.end(), _settings->GetScreenOutput()[idx_field] ) ) ) {

            std::stringstream ss;
            ss << TextProcessingToolbox::DoubleToScientificNotation( _screenOutputFields[idx_field] );
            tmp = ss.str();
            tmp.erase( std::remove( tmp.begin(), tmp.end(), '+' ), tmp.end() );    // removing the '+' sign
        }

        if( strLen > tmp.size() )    // Padding
            tmp.insert( 0, strLen - tmp.size(), paddingChar );
        else if( strLen < tmp.size() )    // Cutting
            tmp.resize( strLen );

        lineToPrint += tmp + " |";
    }
    if( _settings->GetScreenOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetScreenOutputFrequency() == 0 ) {
        log->info( lineToPrint );
    }
    else if( idx_iter == _nIter - 1 ) {    // Always print last iteration
        log->info( lineToPrint );
    }
}

void SNSolverHPC::PrepareHistoryOutput() {
    unsigned n_probes = 4;

    unsigned nFields = (unsigned)_settings->GetNHistoryOutput();

    _historyOutputFieldNames.resize( nFields );
    _historyOutputFields.resize( nFields );

    // Prepare all output Fields ==> Specified in option SCREEN_OUTPUT
    for( unsigned idx_field = 0; idx_field < nFields; idx_field++ ) {
        // Prepare all Output Fields per group

        // Different procedure, depending on the Group...
        switch( _settings->GetHistoryOutput()[idx_field] ) {
            case MASS: _historyOutputFieldNames[idx_field] = "Mass"; break;
            case ITER: _historyOutputFieldNames[idx_field] = "Iter"; break;
            case WALL_TIME: _historyOutputFieldNames[idx_field] = "Wall_time_[s]"; break;
            case RMS_FLUX: _historyOutputFieldNames[idx_field] = "RMS_flux"; break;
            case VTK_OUTPUT: _historyOutputFieldNames[idx_field] = "VTK_out"; break;
            case CSV_OUTPUT: _historyOutputFieldNames[idx_field] = "CSV_out"; break;
            case CUR_OUTFLOW: _historyOutputFieldNames[idx_field] = "Cur_outflow"; break;
            case TOTAL_OUTFLOW: _historyOutputFieldNames[idx_field] = "Total_outflow"; break;
            case CUR_OUTFLOW_P1: _historyOutputFieldNames[idx_field] = "Cur_outflow_P1"; break;
            case TOTAL_OUTFLOW_P1: _historyOutputFieldNames[idx_field] = "Total_outflow_P1"; break;
            case CUR_OUTFLOW_P2: _historyOutputFieldNames[idx_field] = "Cur_outflow_P2"; break;
            case TOTAL_OUTFLOW_P2: _historyOutputFieldNames[idx_field] = "Total_outflow_P2"; break;
            case MAX_OUTFLOW: _historyOutputFieldNames[idx_field] = "Max_outflow"; break;
            case CUR_PARTICLE_ABSORPTION: _historyOutputFieldNames[idx_field] = "Cur_absorption"; break;
            case TOTAL_PARTICLE_ABSORPTION: _historyOutputFieldNames[idx_field] = "Total_absorption"; break;
            case MAX_PARTICLE_ABSORPTION: _historyOutputFieldNames[idx_field] = "Max_absorption"; break;
            case TOTAL_PARTICLE_ABSORPTION_CENTER: _historyOutputFieldNames[idx_field] = "Cumulated_absorption_center"; break;
            case TOTAL_PARTICLE_ABSORPTION_VERTICAL: _historyOutputFieldNames[idx_field] = "Cumulated_absorption_vertical_wall"; break;
            case TOTAL_PARTICLE_ABSORPTION_HORIZONTAL: _historyOutputFieldNames[idx_field] = "Cumulated_absorption_horizontal_wall"; break;
            case PROBE_MOMENT_TIME_TRACE:
                if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) n_probes = 4;
                if( _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) n_probes = 2;
                for( unsigned i = 0; i < n_probes; i++ ) {
                    for( unsigned j = 0; j < 3; j++ ) {
                        _historyOutputFieldNames[idx_field] = "Probe " + std::to_string( i ) + " u_" + std::to_string( j );
                        idx_field++;
                    }
                }
                idx_field--;
                break;
            case VAR_ABSORPTION_GREEN: _historyOutputFieldNames[idx_field] = "Var. absorption green"; break;
            case VAR_ABSORPTION_GREEN_LINE:
                for( unsigned i = 0; i < _settings->GetNumProbingCellsLineHohlraum(); i++ ) {
                    _historyOutputFieldNames[idx_field] = "Probe Green Line " + std::to_string( i );
                    idx_field++;
                }
                idx_field--;
                break;
            default: ErrorMessages::Error( "History output field not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void SNSolverHPC::PrintHistoryOutput( unsigned idx_iter ) {

    auto log = spdlog::get( "tabular" );

    // assemble the line to print
    std::string lineToPrint = "";
    std::string tmp;
    for( int idx_field = 0; idx_field < _settings->GetNHistoryOutput() - 1; idx_field++ ) {
        if( idx_field == 0 ) {
            tmp = std::to_string( _historyOutputFields[idx_field] );    // Iteration count
        }
        else {
            tmp = TextProcessingToolbox::DoubleToScientificNotation( _historyOutputFields[idx_field] );
        }
        lineToPrint += tmp + ",";
    }
    tmp = TextProcessingToolbox::DoubleToScientificNotation( _historyOutputFields[_settings->GetNScreenOutput() - 1] );
    lineToPrint += tmp;    // Last element without comma

    if( _settings->GetHistoryOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetHistoryOutputFrequency() == 0 ) {
        log->info( lineToPrint );
    }
    else if( idx_iter == _nIter - 1 ) {    // Always print last iteration
        log->info( lineToPrint );
    }
}

void SNSolverHPC::DrawPreSolverOutput() {

    // Logger
    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );

    std::string hLine = "--";

    unsigned strLen  = 15;    // max width of one column
    char paddingChar = ' ';

    // Assemble Header for Screen Output
    std::string lineToPrint = "| ";
    std::string tmpLine     = "-----------------";
    for( unsigned idxFields = 0; idxFields < _settings->GetNScreenOutput(); idxFields++ ) {
        std::string tmp = _screenOutputFieldNames[idxFields];

        if( strLen > tmp.size() )    // Padding
            tmp.insert( 0, strLen - tmp.size(), paddingChar );
        else if( strLen < tmp.size() )    // Cutting
            tmp.resize( strLen );

        lineToPrint += tmp + " |";
        hLine += tmpLine;
    }
    log->info( "---------------------------- Solver Starts -----------------------------" );
    log->info( "| The simulation will run for {} iterations.", _nIter );
    log->info( "| The spatial grid contains {} cells.", _nCells );
    if( _settings->GetSolverName() != PN_SOLVER && _settings->GetSolverName() != CSD_PN_SOLVER ) {
        log->info( "| The velocity grid contains {} points.", _nq );
    }
    log->info( hLine );
    log->info( lineToPrint );
    log->info( hLine );

    std::string lineToPrintCSV = "";
    for( int idxFields = 0; idxFields < _settings->GetNHistoryOutput() - 1; idxFields++ ) {
        std::string tmp = _historyOutputFieldNames[idxFields];
        lineToPrintCSV += tmp + ",";
    }
    lineToPrintCSV += _historyOutputFieldNames[_settings->GetNHistoryOutput() - 1];
    logCSV->info( lineToPrintCSV );
}

void SNSolverHPC::DrawPostSolverOutput() {

    // Logger
    auto log = spdlog::get( "event" );

    std::string hLine = "--";

    unsigned strLen  = 10;    // max width of one column
    char paddingChar = ' ';

    // Assemble Header for Screen Output
    std::string lineToPrint = "| ";
    std::string tmpLine     = "------------";
    for( unsigned idxFields = 0; idxFields < _settings->GetNScreenOutput(); idxFields++ ) {
        std::string tmp = _screenOutputFieldNames[idxFields];

        if( strLen > tmp.size() )    // Padding
            tmp.insert( 0, strLen - tmp.size(), paddingChar );
        else if( strLen < tmp.size() )    // Cutting
            tmp.resize( strLen );

        lineToPrint += tmp + " |";
        hLine += tmpLine;
    }
    log->info( hLine );
#ifndef BUILD_TESTING
    log->info( "| The volume output files have been stored at " + _settings->GetOutputFile() );
    log->info( "| The log files have been stored at " + _settings->GetLogDir() + _settings->GetLogFile() );
#endif
    log->info( "--------------------------- Solver Finished ----------------------------" );
}

unsigned SNSolverHPC::Idx2D( unsigned idx1, unsigned idx2, unsigned len2 ) { return idx1 * len2 + idx2; }

unsigned SNSolverHPC::Idx3D( unsigned idx1, unsigned idx2, unsigned idx3, unsigned len2, unsigned len3 ) {
    return ( idx1 * len2 + idx2 ) * len3 + idx3;
}

void SNSolverHPC::WriteVolumeOutput( unsigned idx_iter ) {
    unsigned nGroups = (unsigned)_settings->GetNVolumeOutput();
    if( ( _settings->GetVolumeOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
        ( idx_iter == _nIter - 1 ) /* need sol at last iteration */ ) {
        for( unsigned idx_group = 0; idx_group < nGroups; idx_group++ ) {
            switch( _settings->GetVolumeOutput()[idx_group] ) {
                case MINIMAL:
                    // for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                    _outputFields[idx_group][0] = _scalarFlux;    //[idx_cell];
                    //}
                    break;

                case MOMENTS:
#pragma omp parallel for
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = 0.0;
                        _outputFields[idx_group][1][idx_cell] = 0.0;
                        for( unsigned idx_moments = 0; idx_moments < _nOutputMoments; idx_moments++ ) {
                            _outputFields[idx_group][0][idx_cell] = 0.0;
                            _outputFields[idx_group][1][idx_cell] = 0.0;
                            _outputFields[idx_group][2][idx_cell] = 0.0;
                            // for( unsigned idx_sys = _startSysIdx; idx_sys < _endSysIdx; idx_sys++ ) {    // TODO
                            //     _outputFields[idx_group][0][idx_cell] +=
                            //         _quadPts[Idx2D( idx_sys, 0, _nDim )] * _sol[Idx2D( idx_cell, idx_sys, _localNSys )] * _quadWeights[idx_sys];
                            //     _outputFields[idx_group][1][idx_cell] +=
                            //         _quadPts[Idx2D( idx_sys, 1, _nDim )] * _sol[Idx2D( idx_cell, idx_sys, _localNSys )] * _quadWeights[idx_sys];
                            // }
                        }
                    }
                    break;

                default: ErrorMessages::Error( "Volume Output Group not defined for HPC SN Solver!", CURRENT_FUNCTION ); break;
            }
        }
    }
}

void SNSolverHPC::PrepareVolumeOutput() {
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
                _outputFieldNames[idx_group][0] = "scalar flux";
                break;
            case MOMENTS:
                // As many entries as there are moments in the system
                _outputFields[idx_group].resize( _nOutputMoments );
                _outputFieldNames[idx_group].resize( _nOutputMoments );

                for( unsigned idx_moment = 0; idx_moment < _nOutputMoments; idx_moment++ ) {
                    _outputFieldNames[idx_group][idx_moment] = std::string( "u_" + std::to_string( idx_moment ) );
                    _outputFields[idx_group][idx_moment].resize( _nCells );
                }
                break;

            default: ErrorMessages::Error( "Volume Output Group not defined for HPC SN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

void SNSolverHPC::SetGhostCells() {
    if( _settings->GetProblemName() == PROBLEM_Lattice ) {
        // #pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
            if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NEUMANN || _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) {
                _ghostCells[idx_cell] = std::vector<double>( _localNSys, 0.0 );
            }
        }
    }
    else if( _settings->GetProblemName() == PROBLEM_HalfLattice ) {    // HALF LATTICE NOT WORKING
        ErrorMessages::Error( "Test case does not work with MPI", CURRENT_FUNCTION );
    }
    else if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) {

        auto nodes = _mesh->GetNodes();
        double tol = 1e-12;    // For distance to boundary

        // #pragma omp parallel for
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {

            if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NEUMANN || _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) {
                _ghostCells[idx_cell] = std::vector<double>( _localNSys, 0.0 );

                auto localCellNodes = _mesh->GetCells()[idx_cell];

                for( unsigned idx_node = 0; idx_node < _mesh->GetNumNodesPerCell(); idx_node++ ) {    // Check if corner node is in this cell
                    if( nodes[localCellNodes[idx_node]][0] < -0.65 + tol ) {                          // close to 0 => left boundary
                        for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
                            if( _quadPts[Idx2D( idx_sys, 0, _nDim )] > 0.0 ) _ghostCells[idx_cell][idx_sys] = 1.0;
                        }
                        break;
                    }
                    else if( nodes[localCellNodes[idx_node]][0] > 0.65 - tol ) {    // right boundary
                        for( unsigned idx_sys = 0; idx_sys < _localNSys; idx_sys++ ) {
                            if( _quadPts[Idx2D( idx_sys, 0, _nDim )] < 0.0 ) _ghostCells[idx_cell][idx_sys] = 1.0;
                        }
                        break;
                    }
                    else if( nodes[localCellNodes[idx_node]][1] < -0.65 + tol ) {    // lower boundary
                        break;
                    }
                    else if( nodes[localCellNodes[idx_node]][1] > 0.65 - tol ) {    // upper boundary
                        break;
                    }
                    else if( idx_node == _mesh->GetNumNodesPerCell() - 1 ) {
                        ErrorMessages::Error( " Problem with ghost cell setup and  boundary of this mesh ", CURRENT_FUNCTION );
                    }
                }
            }
        }
    }
    else if( _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) {
        ErrorMessages::Error( "Test case does not work with MPI", CURRENT_FUNCTION );
    }
}

void SNSolverHPC::SetProbingCellsLineGreen() {

    if( _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) {
        double verticalLineWidth   = std::abs( _cornerUpperLeftGreen[1] - _cornerLowerLeftGreen[1] );
        double horizontalLineWidth = std::abs( _cornerUpperLeftGreen[0] - _cornerUpperRightGreen[0] );

        // double dx = 2 * ( horizontalLineWidth + verticalLineWidth ) / ( (double)_nProbingCellsLineGreen );

        unsigned nHorizontalProbingCells =
            (unsigned)std::ceil( _nProbingCellsLineGreen * ( horizontalLineWidth / ( horizontalLineWidth + verticalLineWidth ) ) );
        unsigned nVerticalProbingCells = _nProbingCellsLineGreen - nHorizontalProbingCells;

        _probingCellsLineGreen = std::vector<unsigned>( _nProbingCellsLineGreen );

        // Sample points on each side of the rectangle
        std::vector<unsigned> side3 = linspace2D( _cornerLowerRightGreen, _cornerUpperRightGreen, nVerticalProbingCells );
        std::vector<unsigned> side4 = linspace2D( _cornerUpperRightGreen, _cornerUpperLeftGreen, nHorizontalProbingCells );

        //  Combine the points from each side
        _probingCellsLineGreen.insert( _probingCellsLineGreen.end(), side3.begin(), side3.end() );
        _probingCellsLineGreen.insert( _probingCellsLineGreen.end(), side4.begin(), side4.end() );
    }
    else if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) {
        double verticalLineWidth   = std::abs( _cornerUpperLeftGreen[1] - _cornerLowerLeftGreen[1] );
        double horizontalLineWidth = std::abs( _cornerUpperLeftGreen[0] - _cornerUpperRightGreen[0] );

        // double dx = 2 * ( horizontalLineWidth + verticalLineWidth ) / ( (double)_nProbingCellsLineGreen );

        unsigned nHorizontalProbingCells =
            (unsigned)std::ceil( _nProbingCellsLineGreen / 2 * ( horizontalLineWidth / ( horizontalLineWidth + verticalLineWidth ) ) );
        unsigned nVerticalProbingCells = _nProbingCellsLineGreen - nHorizontalProbingCells;

        _probingCellsLineGreen = std::vector<unsigned>( _nProbingCellsLineGreen );

        std::vector<double> p1 = { _cornerUpperLeftGreen[0] + _thicknessGreen / 2.0, _cornerUpperLeftGreen[1] - _thicknessGreen / 2.0 };
        std::vector<double> p2 = { _cornerLowerLeftGreen[0] + _thicknessGreen / 2.0, _cornerLowerLeftGreen[1] + _thicknessGreen / 2.0 };
        std::vector<double> p3 = { _cornerUpperRightGreen[0] - _thicknessGreen / 2.0, _cornerUpperRightGreen[1] - _thicknessGreen / 2.0 };
        std::vector<double> p4 = { _cornerLowerRightGreen[0] - _thicknessGreen / 2.0, _cornerLowerRightGreen[1] + _thicknessGreen / 2.0 };

        // Sample points on each side of the rectangle
        std::vector<unsigned> side1 = linspace2D( p1, p2, nVerticalProbingCells );
        std::vector<unsigned> side2 = linspace2D( p2, p3, nHorizontalProbingCells );
        std::vector<unsigned> side3 = linspace2D( p3, p4, nVerticalProbingCells );
        std::vector<unsigned> side4 = linspace2D( p4, p1, nHorizontalProbingCells );

        //  Combine the points from each side
        _probingCellsLineGreen.insert( _probingCellsLineGreen.end(), side1.begin(), side1.end() );
        _probingCellsLineGreen.insert( _probingCellsLineGreen.end(), side2.begin(), side2.end() );
        _probingCellsLineGreen.insert( _probingCellsLineGreen.end(), side3.begin(), side3.end() );
        _probingCellsLineGreen.insert( _probingCellsLineGreen.end(), side4.begin(), side4.end() );
    }
}

void SNSolverHPC::ComputeQOIsGreenProbingLine() {
    double verticalLineWidth   = std::abs( _cornerUpperLeftGreen[1] - _cornerLowerLeftGreen[1] - _thicknessGreen );
    double horizontalLineWidth = std::abs( _cornerUpperLeftGreen[0] - _cornerUpperRightGreen[0] - _thicknessGreen );

    double dl    = 2 * ( horizontalLineWidth + verticalLineWidth ) / ( (double)_nProbingCellsLineGreen );
    double area  = dl * _thicknessGreen;
    double a_g   = 0;
    double l_max = _nProbingCellsLineGreen * dl;

    for( unsigned i = 0; i < _nProbingCellsLineGreen; i++ ) {    // Loop over probing cells
        _absorptionValsIntegrated[i] =
            ( _sigmaT[_probingCellsLineGreen[i]] - _sigmaS[_probingCellsLineGreen[i]] ) * _scalarFlux[_probingCellsLineGreen[i]] * area;
        a_g += _absorptionValsIntegrated[i] / (double)_nProbingCellsLineGreen;
    }
    for( unsigned i = 0; i < _nProbingCellsLineGreen; i++ ) {    // Loop over probing cells
        _varAbsorptionValsIntegrated[i] = dl / l_max * ( a_g - _absorptionValsIntegrated[i] ) * ( a_g - _absorptionValsIntegrated[i] );
    }
}

std::vector<unsigned> SNSolverHPC::linspace2D( const std::vector<double>& start, const std::vector<double>& end, unsigned num_points ) {
    /**
     * Generate a 2D linspace based on the start and end points with a specified number of points.
     *
     * @param start vector of starting x and y coordinates
     * @param end vector of ending x and y coordinates
     * @param num_points number of points to generate
     *
     * @return vector of unsigned integers representing the result
     */

    std::vector<unsigned> result;
    result.resize( num_points );
    double stepX = ( end[0] - start[0] ) / ( num_points - 1 );
    double stepY = ( end[1] - start[1] ) / ( num_points - 1 );
#pragma omp parallel for
    for( unsigned i = 0; i < num_points; ++i ) {
        double x = start[0] + i * stepX;
        double y = start[1] + i * stepY;

        result[i] = _mesh->GetCellOfKoordinate( x, y );
    }

    return result;
}

void SNSolverHPC::ComputeCellsPerimeterLattice() {
    double l_1    = 1.5;    // perimeter 1
    double l_2    = 2.5;    // perimeter 2
    auto nodes    = _mesh->GetNodes();
    auto cells    = _mesh->GetCells();
    auto cellMids = _mesh->GetCellMidPoints();
    auto normals  = _mesh->GetNormals();
    auto neigbors = _mesh->GetNeighbours();

    _isPerimeterLatticeCell1.resize( _mesh->GetNumCells(), false );
    _isPerimeterLatticeCell2.resize( _mesh->GetNumCells(), false );

    for( unsigned idx_cell = 0; idx_cell < _mesh->GetNumCells(); ++idx_cell ) {
        if( abs( cellMids[idx_cell][0] ) < l_1 && abs( cellMids[idx_cell][1] ) < l_1 ) {
            // Cell is within perimeter
            for( unsigned idx_nbr = 0; idx_nbr < _mesh->GetNumNodesPerCell(); ++idx_nbr ) {
                if( neigbors[idx_cell][idx_nbr] == _mesh->GetNumCells() ) {
                    continue;    // Skip boundary - ghost cells
                }

                if( abs( ( cellMids[neigbors[idx_cell][idx_nbr]][0] ) > l_1 && abs( cellMids[idx_cell][0] ) < l_1 ) ||
                    abs( ( cellMids[neigbors[idx_cell][idx_nbr]][1] ) > l_1 && abs( cellMids[idx_cell][1] ) < l_1 ) ) {
                    // neighbor is outside perimeter
                    _cellsLatticePerimeter1[idx_cell].push_back( idx_nbr );
                    _isPerimeterLatticeCell1[idx_cell] = true;
                }
            }
        }
        if( abs( cellMids[idx_cell][0] ) < l_2 && abs( cellMids[idx_cell][1] ) < l_2 && abs( cellMids[idx_cell][0] ) > l_1 &&
            abs( cellMids[idx_cell][1] ) > l_1 ) {
            // Cell is within perimeter
            for( unsigned idx_nbr = 0; idx_nbr < _mesh->GetNumNodesPerCell(); ++idx_nbr ) {
                if( neigbors[idx_cell][idx_nbr] == _mesh->GetNumCells() ) {
                    continue;    // Skip boundary - ghost cells
                }
                if( abs( ( cellMids[neigbors[idx_cell][idx_nbr]][0] ) > l_2 && abs( cellMids[idx_cell][0] ) < l_2 ) ||
                    abs( ( cellMids[neigbors[idx_cell][idx_nbr]][1] ) > l_2 && abs( cellMids[idx_cell][1] ) < l_2 ) ) {
                    // neighbor is outside perimeter
                    _cellsLatticePerimeter2[idx_cell].push_back( idx_nbr );
                    _isPerimeterLatticeCell2[idx_cell] = true;
                }
            }
        }
    }
}