#include "solvers/snsolver_hpc.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "kernels/scatteringkernelbase.hpp"
#include "problems/problembase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

SNSolverHPC::SNSolverHPC( Config* settings ) {
    _settings = settings;
    _currTime = 0.0;

    // Create Mesh
    _mesh = LoadSU2MeshFromFile( settings );
    _settings->SetNCells( _mesh->GetNumCells() );
    auto quad = QuadratureBase::Create( settings );
    _settings->SetNQuadPoints( quad->GetNq() );

    _problem = ProblemBase::Create( _settings, _mesh, quad );

    _nCells = _mesh->GetNumCells();
    _nNbr   = _mesh->GetNumNodesPerCell();
    _nDim   = _mesh->GetDim();
    _nNodes = _mesh->GetNumNodes();

    _nq   = quad->GetNq();
    _nSys = _nq;

    _spatialOrder  = _settings->GetSpatialOrder();
    _temporalOrder = _settings->GetTemporalOrder();

    _areas                  = std::vector<double>( _nCells );
    _normals                = std::vector<double>( _nCells * _nNbr * _nDim );
    _neighbors              = std::vector<unsigned>( _nCells * _nNbr );
    _cellMidPoints          = std::vector<double>( _nCells * _nDim );
    _interfaceMidPoints     = std::vector<double>( _nCells * _nNbr * _nDim );
    _cellBoundaryTypes      = std::vector<BOUNDARY_TYPE>( _nCells );
    _relativeInterfaceMidPt = std::vector<double>( _nCells * _nNbr * _nDim );
    // Slope
    _solDx   = std::vector<double>( _nCells * _nSys * _nDim );
    _limiter = std::vector<double>( _nCells * _nSys );

    // Physics
    _sigmaS           = std::vector<double>( _nCells );
    _sigmaT           = std::vector<double>( _nCells );
    _source           = std::vector<double>( _nCells * _nSys );
    _scatteringKernel = std::vector<double>( _nSys * _nSys );

    // Quadrature
    _quadPts               = std::vector<double>( _nSys * _nDim );
    _quadWeights           = std::vector<double>( _nSys );
    _scatteringKernel      = std::vector<double>( _nSys * _nSys );
    _quadratureYReflection = std::vector<unsigned>( _nSys );

    // Solution
    _sol    = std::vector<double>( _nCells * _nSys );
    _solNew = std::vector<double>( _nCells * _nSys );

    _scalarFlux              = std::vector<double>( _nCells );
    _scalarFluxNew           = std::vector<double>( _nCells );
    _localMaxOrdinateOutflow = std::vector<double>( _nCells );

    auto areas           = _mesh->GetCellAreas();
    auto neighbors       = _mesh->GetNeighbours();
    auto normals         = _mesh->GetNormals();
    auto nodes           = _mesh->GetNodes();
    auto cellNodes       = _mesh->GetCells();
    auto cellMidPts      = _mesh->GetCellMidPoints();
    auto interfaceMidPts = _mesh->GetInterfaceMidPoints();
    auto boundaryTypes   = _mesh->GetBoundaryTypes();

    for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
        _areas[idx_cell] = areas[idx_cell];
    }

    _dT    = ComputeTimeStep( _settings->GetCFL() );
    _nIter = unsigned( _settings->GetTEnd() / _dT ) + 1;

    auto quadPoints  = quad->GetPoints();
    auto quadWeights = quad->GetWeights();

    auto initialCondition = _problem->SetupIC();
    auto sigmaT           = _problem->GetTotalXS( Vector( _nIter, 0.0 ) );
    auto sigmaS           = _problem->GetScatteringXS( Vector( _nIter, 0.0 ) );
    auto source           = _problem->GetExternalSource( Vector( _nIter, 0.0 ) );

    SetGhostCells();    // To be changed

// Copy to everything to solver
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
            }
        }

        _sigmaS[idx_cell]     = sigmaS[0][idx_cell];
        _sigmaT[idx_cell]     = sigmaT[0][idx_cell];
        _scalarFlux[idx_cell] = 0;
        for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {
            _source[Idx2D( idx_cell, idx_sys, _nSys )] = source[0][idx_cell][0];    // CAREFUL HERE hardcoded to isotropic source
            _sol[Idx2D( idx_cell, idx_sys, _nSys )]    = initialCondition[idx_cell][idx_sys];
            _solNew[Idx2D( idx_cell, idx_sys, _nSys )] = initialCondition[idx_cell][idx_sys];
            _scalarFlux[idx_cell] += _solNew[Idx2D( idx_cell, idx_sys, _nSys )] * _quadWeights[idx_sys];
        }
    }

    ScatteringKernel* k   = ScatteringKernel::CreateScatteringKernel( settings->GetKernelName(), quad );
    auto scatteringKernel = k->GetScatteringKernel();
    for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {
        for( unsigned idx_dim = 0; idx_dim < _nDim; idx_dim++ ) {
            _quadPts[Idx2D( idx_sys, idx_dim, _nDim )] = quadPoints[idx_sys][idx_dim];
        }
        _quadWeights[idx_sys] = quadWeights[idx_sys];

        for( unsigned idx_sys2 = 0; idx_sys2 < _nSys; idx_sys2++ ) {
            _scatteringKernel[Idx2D( idx_sys, idx_sys2, _nSys )] = scatteringKernel( idx_sys, idx_sys2 );
        }
    }

    PrepareScreenOutput();     // Screen Output
    PrepareHistoryOutput();    // History Output

    delete quad;
    delete k;

    // Initialiye QOIS
    _mass                    = 0;
    _rmsFlux                 = 0;
    _curAbsorptionLattice    = 0;
    _totalAbsorptionLattice  = 0;
    _curMaxAbsorptionLattice = 0;
}

void SNSolverHPC::Solve() {

    // --- Preprocessing ---
    PrepareVolumeOutput();

    DrawPreSolverOutput();

    // Create Backup solution for Runge Kutta
    // std::vector<double> solRK0 = _sol;

    auto start = std::chrono::high_resolution_clock::now();    // Start timing

    std::chrono::duration<double> duration;
    // Loop over energies (pseudo-time of continuous slowing down approach)
    for( unsigned iter = 0; iter < _nIter; iter++ ) {
        if( iter == _nIter - 1 ) {    // last iteration
            _dT = _settings->GetTEnd() - iter * _dT;
        }
        // if( _temporalOrder == 2 ) solRK0 = _sol;

        // for( unsigned rkStep = 0; rkStep < _temporalOrder; ++rkStep ) {
        //  --- Prepare Boundaries and temp variables

        ( _spatialOrder == 2 ) ? FVMUpdateOrder2() : FVMUpdateOrder1();

        // --- Update Solution within Runge Kutta Stages

        //}

        // --- Wall time measurement
        duration  = std::chrono::high_resolution_clock::now() - start;
        _currTime = std::chrono::duration_cast<std::chrono::duration<double>>( duration ).count();

        // --- Runge Kutta Timestep ---
        // if( _temporalOrder == 2 ) RKUpdate( solRK0, _sol );

        // --- Write Output ---
        WriteVolumeOutput( iter );
        WriteScalarOutput( iter );

        // --- Update Scalar Fluxes

        // --- Print Output ---
        PrintScreenOutput( iter );
        PrintHistoryOutput( iter );
        PrintVolumeOutput( iter );
    }

    // --- Postprocessing ---

    DrawPostSolverOutput();
}

void SNSolverHPC::RKUpdate( std::vector<double>& sol_0, std::vector<double>& sol_rk ) {
#pragma omp parallel for
    for( unsigned i = 0; i < _nCells * _nSys; ++i ) {
        _sol[i] = 0.5 * ( sol_0[i] + sol_rk[i] );
    }
}

void SNSolverHPC::PrintVolumeOutput() const { ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh ); }

void SNSolverHPC::PrintVolumeOutput( int idx_iter ) const {
    if( _settings->GetVolumeOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) {
        ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( idx_iter ), _outputFields, _outputFieldNames, _mesh );    // slow
    }
    if( idx_iter == (int)_nIter - 1 ) {    // Last iteration write without suffix.
        ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh );
    }
}

void SNSolverHPC::FVMUpdateOrder2() {
    double const eps = 1e-10;
    _mass            = 0.0;
    _rmsFlux         = 0.0;
#pragma omp parallel
    {
        double localMass    = 0.0;
        double localRMSFlux = 0.0;

        double localMaxAbsorptionLattice  = _curMaxAbsorptionLattice;
        double localCurrAbsorptionLattice = 0.0;

        double localMaxOrdinatewiseOutflow = _curMaxOrdinateOutflow;
        double localCurrScalarOutflow      = 0.0;

        // Loop over all spatial cells

        double solL;
        double solR;
        double inner;
        double minSol;
        double maxSol;
        double r;

        unsigned idx_nbr_glob = 0;

        double solInterfaceAvg = 0.0;

        double solDx;
        double solDy;
        double limiter;
        std::vector<double> gaussPoints( _nNbr );

#pragma omp for nowait
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {

            // Reset temporary variable
            for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {
                _solNew[Idx2D( idx_cell, idx_sys, _nSys )] = 0.0;
            }

            // Fluxes

            // First order at boundary
            if( _cellBoundaryTypes[idx_cell] != BOUNDARY_TYPE::NONE ) {
                // Fluxes
                for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {
                    if( _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )] == _nCells ) {
                        for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {
                            inner = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                    _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];

                            if( inner > 0 ) {
                                _solNew[Idx2D( idx_cell, idx_sys, _nSys )] += inner * _sol[Idx2D( idx_cell, idx_sys, _nSys )];
                            }
                            else {

                                double ghostCellValue = ( _ghostCellsReflectingY[idx_cell] )
                                                            ? _sol[Idx2D( idx_cell, _quadratureYReflection[idx_sys], _nSys )]    // Relecting boundary
                                                            : _ghostCells[idx_cell][idx_sys];                                    // fixed boundary

                                _solNew[Idx2D( idx_cell, idx_sys, _nSys )] += inner * ghostCellValue;
                            }
                        }
                    }
                    else {
                        idx_nbr_glob = _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )];    // global idx of neighbor cell
                        for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {
                            inner = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                    _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];
                            _solNew[Idx2D( idx_cell, idx_sys, _nSys )] +=
                                ( inner > 0 ) ? inner * _sol[Idx2D( idx_cell, idx_sys, _nSys )] : inner * _sol[Idx2D( idx_nbr_glob, idx_sys, _nSys )];
                        }
                    }
                }
            }

            else {    // Second order
                for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {

                    // Slope computation
                    solDx   = 0.0;
                    solDy   = 0.0;
                    limiter = 1.0;

                    minSol = _sol[Idx2D( idx_cell, idx_sys, _nSys )];
                    maxSol = minSol;
                    for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {
                        //  compute derivative by summing over cell boundary using Green Gauss Theorem
                        idx_nbr_glob    = _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )];
                        solInterfaceAvg = 0.5 * ( _sol[Idx2D( idx_cell, idx_sys, _nSys )] + _sol[Idx2D( idx_nbr_glob, idx_sys, _nSys )] );

                        solDx += solInterfaceAvg * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )];
                        solDy += solInterfaceAvg * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];

                        // Compute ptswise max and minimum solultion values of current and neighbor cells
                        maxSol = std::max( _sol[Idx2D( idx_nbr_glob, idx_sys, _nSys )], maxSol );
                        minSol = std::min( _sol[Idx2D( idx_nbr_glob, idx_sys, _nSys )], minSol );
                    }
                    solDx /= _areas[idx_cell];
                    solDy /= _areas[idx_cell];

                    for( unsigned idx_nbr = 0; idx_nbr < _nNbr; idx_nbr++ ) {

                        // Compute value at interface midpoint, called gaussPt
                        gaussPoints[idx_nbr] = solDx * _relativeInterfaceMidPt[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                               solDy * _relativeInterfaceMidPt[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];

                        // BARTH-JESPERSEN LIMITER
                        // r = ( gaussPoints[idx_nbr] > 0 ) ? ( maxSol - _sol[Idx2D( idx_cell, idx_sys, _nSys )] ) / ( gaussPoints[idx_nbr] + eps )
                        //                                 : ( minSol - _sol[Idx2D( idx_cell, idx_sys, _nSys )] ) / ( gaussPoints[idx_nbr] - eps );
                        //                        limiter = std::min( std::min( r, 1.0 ), limiter );

                        // VENKATAKRISHNAN LIMITER
                        double delta1Max = maxSol - _sol[Idx2D( idx_cell, idx_sys, _nSys )];
                        double delta1Min = minSol - _sol[Idx2D( idx_cell, idx_sys, _nSys )];

                        r = ( gaussPoints[idx_nbr] > 0 ) ? ( ( delta1Max * delta1Max + _areas[idx_cell] ) * gaussPoints[idx_nbr] +
                                                             2 * gaussPoints[idx_nbr] * gaussPoints[idx_nbr] * delta1Max ) /
                                                               ( delta1Max * delta1Max + 2 * gaussPoints[idx_nbr] * gaussPoints[idx_nbr] +
                                                                 delta1Max * gaussPoints[idx_nbr] + _areas[idx_cell] ) /
                                                               ( gaussPoints[idx_nbr] + eps )
                                                         : ( ( delta1Min * delta1Min + _areas[idx_cell] ) * gaussPoints[idx_nbr] +
                                                             2 * gaussPoints[idx_nbr] * gaussPoints[idx_nbr] * delta1Min ) /
                                                               ( delta1Min * delta1Min + 2 * gaussPoints[idx_nbr] * gaussPoints[idx_nbr] +
                                                                 delta1Min * gaussPoints[idx_nbr] + _areas[idx_cell] ) /
                                                               ( gaussPoints[idx_nbr] - eps );

                        r       = ( std::abs( gaussPoints[idx_nbr] ) < eps ) ? 1 : r;
                        limiter = std::min( r, limiter );
                    }

                    for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {
                        // store flux contribution on psiNew_sigmaS to save memory
                        idx_nbr_glob = _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )];    // global idx of neighbor cell
                        inner        = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];

                        _solNew[Idx2D( idx_cell, idx_sys, _nSys )] +=
                            ( inner > 0 ) ? inner * _sol[Idx2D( idx_cell, idx_sys, _nSys )] + limiter * gaussPoints[idx_nbr]
                                          : inner * _sol[Idx2D( idx_nbr_glob, idx_sys, _nSys )] + limiter * gaussPoints[idx_nbr];
                    }
                }
            }
            // Upate
            for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {
                _solNew[Idx2D( idx_cell, idx_sys, _nSys )] = _sol[Idx2D( idx_cell, idx_sys, _nSys )] -
                                                             ( _dT / _areas[idx_cell] ) * _solNew[Idx2D( idx_cell, idx_sys, _nSys )] -
                                                             _dT * _sigmaT[idx_cell] * _sol[Idx2D( idx_cell, idx_sys, _nSys )];

                // In-scattering Term

                _solNew[Idx2D( idx_cell, idx_sys, _nSys )] +=
                    _dT * _sigmaS[idx_cell] * _scalarFlux[idx_cell] / ( 4 * M_PI );    // Isotropic scattering

                // Source Term
                _solNew[Idx2D( idx_cell, idx_sys, _nSys )] += _dT * _source[Idx2D( idx_cell, idx_sys, _nSys )];
            }
            // --- Iter Postprocessing ---

            _scalarFluxNew[idx_cell] = 0;
            for( unsigned idx_sys = 0; idx_sys < _nSys; ++idx_sys ) {
                _scalarFluxNew[idx_cell] += _solNew[Idx2D( idx_cell, idx_sys, _nSys )] * _quadWeights[idx_sys];
                _sol[Idx2D( idx_cell, idx_sys, _nSys )] = _solNew[Idx2D( idx_cell, idx_sys, _nSys )];    // reset sol
            }

            localMass += _scalarFluxNew[idx_cell] * _areas[idx_cell];
            localRMSFlux += ( _scalarFluxNew[idx_cell] - _scalarFlux[idx_cell] ) * ( _scalarFluxNew[idx_cell] - _scalarFlux[idx_cell] );
            _scalarFlux[idx_cell] = _scalarFluxNew[idx_cell];    // reset flux

            //_problem->ComputeCurrentAbsorptionLattice( _scalarFlux );
            if( IsAbsorptionLattice( _cellMidPoints[Idx2D( idx_cell, 0, _nDim )], _cellMidPoints[Idx2D( idx_cell, 1, _nDim )] ) ) {
                localCurrAbsorptionLattice += _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * _areas[idx_cell];

                if( localMaxAbsorptionLattice < _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) )

                    localMaxAbsorptionLattice = _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] );    // Need to be a private area
            }

            if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NEUMANN && !_ghostCellsReflectingY[idx_cell] ) {
                // Iterate over face cell faces
                double currOrdinatewiseOutflow = 0.0;

                for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {
                    // Find face that points outward
                    if( _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )] == _nCells ) {
                        // Iterate over transport directions
                        for( unsigned idx_sys = 0; idx_sys < _nSys; ++idx_sys ) {
                            inner = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                    _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];
                            // Find outward facing transport directions
                            if( inner > 0.0 ) {
                                localCurrScalarOutflow +=
                                    inner * _sol[Idx2D( idx_cell, idx_sys, _nSys )] * _quadWeights[idx_sys];    // Integrate flux

                                currOrdinatewiseOutflow = _sol[Idx2D( idx_cell, idx_sys, _nSys )] * inner /
                                                          sqrt( ( _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] *
                                                                      _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                                                  _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )] *
                                                                      _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )] ) );

                                if( currOrdinatewiseOutflow > localMaxOrdinatewiseOutflow ) {
                                    localMaxOrdinatewiseOutflow = currOrdinatewiseOutflow;
                                }
                            }
                        }
                    }
                }
            }
        }
#pragma omp critical
        {
            if( localMaxOrdinatewiseOutflow > _curMaxOrdinateOutflow ) {
                _curMaxOrdinateOutflow = localMaxOrdinatewiseOutflow;
            }
            // reduction( + : _mass, _rmsFlux )
            _rmsFlux += localRMSFlux;
            _mass += localMass;
            _curAbsorptionLattice += localCurrAbsorptionLattice;
            _curScalarOutflow += localCurrScalarOutflow;

            if( localMaxAbsorptionLattice > _curMaxAbsorptionLattice ) {
                _curMaxAbsorptionLattice = localMaxAbsorptionLattice;
            }
        }
    }
    _rmsFlux = sqrt( _rmsFlux );

    _totalScalarOutflow += _curScalarOutflow * _dT;
    _totalAbsorptionLattice += _curAbsorptionLattice * _dT;
}

void SNSolverHPC::FVMUpdateOrder1() {
    _mass    = 0.0;
    _rmsFlux = 0.0;
#pragma omp parallel
    {
        double localMass    = 0.0;
        double localRMSFlux = 0.0;

        double localMaxAbsorptionLattice  = _curMaxAbsorptionLattice;
        double localCurrAbsorptionLattice = 0.0;

        double localMaxOrdinatewiseOutflow = _curMaxOrdinateOutflow;
        double localCurrScalarOutflow      = 0.0;

        // Loop over all spatial cells

#pragma omp for nowait
        for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
            double solL;
            double solR;
            double inner;
            unsigned nbr_glob;
            // Reset temporary variable
            for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {
                _solNew[Idx2D( idx_cell, idx_sys, _nSys )] = 0.0;
            }

            // Fluxes
            for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {
                if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NEUMANN && _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )] == _nCells ) {
                    for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {
                        inner = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];

                        if( inner > 0 ) {
                            _solNew[Idx2D( idx_cell, idx_sys, _nSys )] += inner * _sol[Idx2D( idx_cell, idx_sys, _nSys )];
                        }
                        else {

                            double ghostCellValue = ( _ghostCellsReflectingY[idx_cell] )
                                                        ? _sol[Idx2D( idx_cell, _quadratureYReflection[idx_sys], _nSys )]    // Relecting boundary
                                                        : _ghostCells[idx_cell][idx_sys];                                    // fixed boundary

                            _solNew[Idx2D( idx_cell, idx_sys, _nSys )] += inner * ghostCellValue;
                        }
                    }
                }
                else {
                    nbr_glob = _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )];    // global idx of neighbor cell
                    for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {
                        inner = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];
                        _solNew[Idx2D( idx_cell, idx_sys, _nSys )] +=
                            ( inner > 0 ) ? inner * _sol[Idx2D( idx_cell, idx_sys, _nSys )] : inner * _sol[Idx2D( nbr_glob, idx_sys, _nSys )];
                    }
                }
            }
            // Upate
            for( unsigned idx_sys = 0; idx_sys < _nSys; idx_sys++ ) {
                _solNew[Idx2D( idx_cell, idx_sys, _nSys )] = _sol[Idx2D( idx_cell, idx_sys, _nSys )] -
                                                             ( _dT / _areas[idx_cell] ) * _solNew[Idx2D( idx_cell, idx_sys, _nSys )] -
                                                             _dT * _sigmaT[idx_cell] * _sol[Idx2D( idx_cell, idx_sys, _nSys )];

                // In-scattering Term

                _solNew[Idx2D( idx_cell, idx_sys, _nSys )] +=
                    _dT * _sigmaS[idx_cell] * _scalarFlux[idx_cell] / ( 4 * M_PI );    // Isotropic scattering

                // Source Term
                _solNew[Idx2D( idx_cell, idx_sys, _nSys )] += _dT * _source[Idx2D( idx_cell, idx_sys, _nSys )];
            }

            // --- Iter Postprocessing ---

            _scalarFluxNew[idx_cell] = 0;
            for( unsigned idx_sys = 0; idx_sys < _nSys; ++idx_sys ) {
                _scalarFluxNew[idx_cell] += _solNew[Idx2D( idx_cell, idx_sys, _nSys )] * _quadWeights[idx_sys];
                _sol[Idx2D( idx_cell, idx_sys, _nSys )] = _solNew[Idx2D( idx_cell, idx_sys, _nSys )];    // reset sol
            }

            localMass += _scalarFluxNew[idx_cell] * _areas[idx_cell];
            localRMSFlux += ( _scalarFluxNew[idx_cell] - _scalarFlux[idx_cell] ) * ( _scalarFluxNew[idx_cell] - _scalarFlux[idx_cell] );
            _scalarFlux[idx_cell] = _scalarFluxNew[idx_cell];    // reset flux

            if( _settings->GetProblemName() == PROBLEM_Lattice || _settings->GetProblemName() == PROBLEM_HalfLattice ) {
                //_problem->ComputeCurrentAbsorptionLattice( _scalarFlux );
                if( IsAbsorptionLattice( _cellMidPoints[Idx2D( idx_cell, 0, _nDim )], _cellMidPoints[Idx2D( idx_cell, 1, _nDim )] ) ) {
                    localCurrAbsorptionLattice += _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * _areas[idx_cell];

                    if( localMaxAbsorptionLattice < _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) )

                        localMaxAbsorptionLattice = _scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] );    // Need to be a private area
                }
            }

            if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NEUMANN ) {
                // Iterate over face cell faces
                double currOrdinatewiseOutflow = 0.0;

                for( unsigned idx_nbr = 0; idx_nbr < _nNbr; ++idx_nbr ) {
                    // Find face that points outward
                    if( _neighbors[Idx2D( idx_cell, idx_nbr, _nNbr )] == _nCells ) {
                        // Iterate over transport directions
                        for( unsigned idx_sys = 0; idx_sys < _nSys; ++idx_sys ) {
                            inner = _quadPts[Idx2D( idx_sys, 0, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                    _quadPts[Idx2D( idx_sys, 1, _nDim )] * _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )];
                            // Find outward facing transport directions
                            if( inner > 0.0 ) {
                                localCurrScalarOutflow +=
                                    inner * _sol[Idx2D( idx_cell, idx_sys, _nSys )] * _quadWeights[idx_sys];    // Integrate flux

                                currOrdinatewiseOutflow = _sol[Idx2D( idx_cell, idx_sys, _nSys )] * inner /
                                                          sqrt( ( _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] *
                                                                      _normals[Idx3D( idx_cell, idx_nbr, 0, _nNbr, _nDim )] +
                                                                  _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )] *
                                                                      _normals[Idx3D( idx_cell, idx_nbr, 1, _nNbr, _nDim )] ) );

                                if( currOrdinatewiseOutflow > localMaxOrdinatewiseOutflow ) {
                                    localMaxOrdinatewiseOutflow = currOrdinatewiseOutflow;
                                }
                            }
                        }
                    }
                }
            }
        }
#pragma omp critical
        {
            if( localMaxOrdinatewiseOutflow > _curMaxOrdinateOutflow ) {
                _curMaxOrdinateOutflow = localMaxOrdinatewiseOutflow;
            }
            // reduction( + : _mass, _rmsFlux )
            _rmsFlux += localRMSFlux;
            _mass += localMass;
            _curAbsorptionLattice += localCurrAbsorptionLattice;
            _curScalarOutflow += localCurrScalarOutflow;

            if( localMaxAbsorptionLattice > _curMaxAbsorptionLattice ) {
                _curMaxAbsorptionLattice = localMaxAbsorptionLattice;
            }
        }
    }
    _rmsFlux = sqrt( _rmsFlux );

    _totalScalarOutflow += _curScalarOutflow * _dT;
    _totalAbsorptionLattice += _curAbsorptionLattice * _dT;
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
    for( unsigned j = 0; j < _nCells; j++ ) {
        double currCharSize = sqrt( _areas[j] );
        if( currCharSize < charSize ) {
            charSize = currCharSize;
        }
    }
    auto log         = spdlog::get( "event" );
    std::string line = "| Smallest characteristic length of a grid cell in this mesh: " + std::to_string( charSize );
    log->info( line );
    line = "| Corresponding maximal time-step: " + std::to_string( cfl * charSize );
    log->info( line );
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
            case MAX_OUTFLOW: _screenOutputFields[idx_field] = _curMaxOrdinateOutflow; break;
            case CUR_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _curAbsorptionLattice; break;
            case TOTAL_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _totalAbsorptionLattice; break;
            case MAX_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _curMaxAbsorptionLattice; break;
            case TOTAL_PARTICLE_ABSORPTION_CENTER: _screenOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumCenter(); break;
            case TOTAL_PARTICLE_ABSORPTION_VERTICAL: _screenOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumVertical(); break;
            case TOTAL_PARTICLE_ABSORPTION_HORIZONTAL: _screenOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumHorizontal(); break;
            case PROBE_MOMENT_TIME_TRACE:
                if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) n_probes = 4;
                if( _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) n_probes = 2;
                for( unsigned i = 0; i < n_probes; i++ ) {
                    _screenOutputFields[idx_field] = probingMoments[i][0];
                    idx_field++;
                }
                idx_field--;
                break;
            case VAR_ABSORPTION_GREEN: _screenOutputFields[idx_field] = _problem->GetVarAbsorptionHohlraumGreen(); break;
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
            case MAX_OUTFLOW: _historyOutputFields[idx_field] = _curMaxOrdinateOutflow; break;
            case CUR_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _curAbsorptionLattice; break;
            case TOTAL_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _totalAbsorptionLattice; break;
            case MAX_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _curMaxAbsorptionLattice; break;
            case TOTAL_PARTICLE_ABSORPTION_CENTER: _historyOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumCenter(); break;
            case TOTAL_PARTICLE_ABSORPTION_VERTICAL: _historyOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumVertical(); break;
            case TOTAL_PARTICLE_ABSORPTION_HORIZONTAL: _historyOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumHorizontal(); break;
            case PROBE_MOMENT_TIME_TRACE:
                if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) n_probes = 4;
                if( _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) n_probes = 2;
                for( unsigned i = 0; i < n_probes; i++ ) {
                    for( unsigned j = 0; j < 3; j++ ) {
                        _historyOutputFields[idx_field] = probingMoments[i][j];
                        idx_field++;
                    }
                }
                idx_field--;
                break;
            case VAR_ABSORPTION_GREEN: _historyOutputFields[idx_field] = _problem->GetVarAbsorptionHohlraumGreen(); break;
            case VAR_ABSORPTION_GREEN_LINE:
                for( unsigned i = 0; i < _settings->GetNumProbingCellsLineHohlraum(); i++ ) {
                    _historyOutputFieldNames[idx_field] = _problem->GetCurrentVarProbeValuesGreenLine()[i];
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
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        _outputFields[idx_group][0][idx_cell] = _scalarFluxNew[idx_cell];
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

            default: ErrorMessages::Error( "Volume Output Group not defined for HPC SN Solver!", CURRENT_FUNCTION ); break;
        }
    }
}

void SNSolverHPC::SetGhostCells() {
    if( _settings->GetProblemName() == PROBLEM_Lattice ) {
        for( unsigned idx_cell = 0; idx_cell < _nCells; idx_cell++ ) {
            if( _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::NEUMANN || _cellBoundaryTypes[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) {
                _ghostCells.insert( { idx_cell, std::vector<double>( _nSys ) } );
                _ghostCellsReflectingY[idx_cell] = false;
            }
        }
    }
    else {    // HALF LATTICE

        double tol = 1e-12;    // For distance to boundary

        QuadratureBase* quad = QuadratureBase::Create( _settings );
        VectorVector vq      = quad->GetPoints();
        unsigned nq          = quad->GetNq();

        if( _settings->GetQuadName() != QUAD_GaussLegendreTensorized2D ) {
            ErrorMessages::Error( "This simplified test case only works with symmetric quadrature orders. Use QUAD_GAUSS_LEGENDRE_TENSORIZED_2D",
                                  CURRENT_FUNCTION );
        }
        {    // Create the symmetry maps for the quadratures

            for( unsigned idx_q = 0; idx_q < nq; idx_q++ ) {
                for( unsigned idx_q2 = 0; idx_q2 < nq; idx_q2++ ) {
                    if( abs( vq[idx_q][0] + vq[idx_q2][0] ) + abs( vq[idx_q][1] - vq[idx_q2][1] ) < tol ) {
                        _quadratureYReflection[idx_q] = idx_q2;
                        break;
                    }
                }
            }
        }

        if( _quadratureYReflection.size() != nq ) {
            ErrorMessages::Error( "Problem with Y symmetry of quadrature of this mesh", CURRENT_FUNCTION );
        }

        auto nodes = _mesh->GetNodes();

        for( unsigned idx_cell = 0; idx_cell < _mesh->GetNumCells(); idx_cell++ ) {
            if( _boundaryCells[idx_cell] == BOUNDARY_TYPE::NEUMANN || _boundaryCells[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) {
                double x = _mesh->GetCellMidPoints()[idx_cell][0];
                double y = _mesh->GetCellMidPoints()[idx_cell][1];

                auto localCellNodes = _mesh->GetCells()[idx_cell];

                _ghostCellsReflectingY[idx_cell] = false;
                for( unsigned idx_node = 0; idx_node < _mesh->GetNumNodesPerCell(); idx_node++ ) {    // Check if corner node is in this cell
                    if( abs( nodes[localCellNodes[idx_node]][1] ) < -3.5 + tol || abs( nodes[localCellNodes[idx_node]][1] ) > 3.5 - tol ) {
                        // upper and lower boundary
                        _ghostCells.insert( { idx_cell, std::vector<double>( nq ) } );
                        break;
                    }
                    else if( abs( nodes[localCellNodes[idx_node]][0] ) < tol ) {    // close to 0 => left boundary
                        _ghostCellsReflectingY[idx_cell] = true;
                        _ghostCells.insert( { idx_cell, std::vector<double>( nq ) } );
                        break;
                    }
                    else {    // right boundary
                        _ghostCells.insert( { idx_cell, std::vector<double>( nq ) } );
                        break;
                    }
                }
            }
        }

        delete quad;
    }
}
