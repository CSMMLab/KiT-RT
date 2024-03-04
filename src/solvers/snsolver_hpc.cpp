#include "solvers/snsolver_hpc.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "kernels/scatteringkernelbase.hpp"
#include "problems/problembase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "quadratures/quadraturebase.hpp"

SNSolverHPC::SNSolverHPC( Config* settings ) {
    _settings = settings;
    _currTime = 0.0;


    // Create Mesh
    auto mesh = LoadSU2MeshFromFile( settings );    
    _settings->SetNCells( mesh->GetNumCells() );
    auto quad = QuadratureBase::Create( settings );
    _settings->SetNQuadPoints( _nq );
    auto problem = ProblemBase::Create( _settings, mesh, quad );

    _nCells     = mesh->GetNumCells();
    _nEdgesPerCell = mesh->GetNumNodesPerCell();
    _spatialDim = mesh->GetDim();   
    _nq         = quad->GetNq();
    _nSystem = _nq; 

    _areas = std::vector<double>( _nCells );
    _normals = std::vector<double>( _nCells *_nEdgesPerCell*_spatialDim );
    _neighbors = std::vector<unsigned>( _nCells *_nEdgesPerCell );
    _cellMidPoints = std::vector<double>( _nCells *_nEdgesPerCell );
    _interfaceMidPoints = std::vector<double>( _nCells *_nEdgesPerCell*_spatialDim );
    
    // Slope
    _solDx    = std::vector<double>( _nCells *_nSystem*_spatialDim );
    _limiter    = std::vector<double>( _nCells *_spatialDim );
    
    // Physics
    _sigmaS   = std::vector<double>( _nCells );
    _sigmaT   = std::vector<double>( _nCells );
    _source   = std::vector<double>( _nCells*_nSystem );
    _scatteringKernel   = std::vector<double>( _nSystem * _nSystem );

    //Quadrature
    _quadPts   = std::vector<double>( _nSystem * _spatialDim );
    _quadWeights   = std::vector<double>( _nSystem );
    _scatteringKernel =std::vector<double>( _nSystem *_nSystem);

    // Solution 
    _sol = std::vector<double>( _nCells*_nSystem );
    _solNew = std::vector<double>( _nCells*_nSystem );
    
    auto areas           = mesh->GetCellAreas();
    auto neighbors       = mesh->GetNeighbours();
    auto normals         = mesh->GetNormals();
    auto cellMidPts      = mesh->GetCellMidPoints();
    auto interfaceMidPts = mesh->GetInterfaceMidPoints();

    _dT = ComputeTimeStep( _settings->GetCFL() );
    _nIter =  ceil( _settings->GetTEnd() / _dT );    // redundancy with nIter (<-Better) ?


    auto quadPoints = quad->GetPoints();
    auto quadWeights    = quad->GetWeights();

    auto initialCondition = problem->SetupIC();
    auto sigmaT = problem->GetTotalXS( Vector(_nIter,0.0) );
    auto sigmaS = problem->GetScatteringXS( Vector(_nIter,0.0) );
    auto source      = problem->GetExternalSource( Vector(_nIter,0.0) );

    // Copy to everything to solver
    for (unsigned idx_cell=0;idx_cell<_nCells; idx_cell ++){
        _areas[idx_cell] =areas[idx_cell];
        for (unsigned idx_nbr=0;idx_nbr<_nEdgesPerCell; idx_nbr ++){
            _neighbors[idx_cell*_nEdgesPerCell + idx_nbr] = neighbors[idx_cell][idx_nbr];
            _cellMidPoints[idx_cell*_nEdgesPerCell + idx_nbr] = cellMidPts[idx_cell][ idx_nbr];

            for (unsigned idx_dim=0;idx_dim<_spatialDim; idx_dim ++){
                _normals[idx_cell*_nEdgesPerCell*_spatialDim + idx_nbr*_spatialDim +idx_dim] =normals[idx_cell][idx_nbr][idx_dim];
                _interfaceMidPoints[idx_cell*_nEdgesPerCell*_spatialDim + idx_nbr*_spatialDim +idx_dim]  =interfaceMidPts[idx_cell][idx_nbr][idx_dim];
            }
        }
        _sigmaS[idx_cell] =sigmaS[0][idx_cell];
        _sigmaT[idx_cell] =sigmaT[0][idx_cell];
        for (unsigned idx_sys=0;idx_sys<_nSystem; idx_sys ++){
            _source[idx_cell*_nSystem + idx_sys] =source[0][idx_cell][0]; //CAREFUL HERE hardcoded to isotropic source
                 _sol[idx_cell*_nSystem + idx_sys] =   initialCondition  [idx_cell][idx_sys];
                 _solNew[idx_cell*_nSystem + idx_sys] =   initialCondition  [idx_cell][idx_sys];
        }

    }   

    ScatteringKernel* k = ScatteringKernel::CreateScatteringKernel( settings->GetKernelName(), quad );
    auto scatteringKernel   = k->GetScatteringKernel();
    for (unsigned idx_sys=0;idx_sys<_nSystem; idx_sys ++){
        for (unsigned idx_sys2=0;idx_sys2<_nSystem; idx_sys2 ++){
            _scatteringKernel[idx_sys*_nSystem +idx_sys2 ] =  scatteringKernel(idx_sys,idx_sys2 );
        }
    }

    PrepareScreenOutput();     // Screen Output
    PrepareHistoryOutput();    // History Output

    delete mesh;
    delete quad;
    delete problem;
    delete k;
}


void SNSolverHPC::Solve() {

    // --- Preprocessing ---
    PrepareVolumeOutput();

    DrawPreSolverOutput();

    // Preprocessing before first pseudo time step
    SolverPreprocessing();
    unsigned rkStages = _settings->GetRKStages();
    // Create Backup solution for Runge Kutta
    VectorVector solRK0 = _sol;

    auto start = std::chrono::high_resolution_clock::now();    // Start timing
    std::chrono::duration<double> duration;
    // Loop over energies (pseudo-time of continuous slowing down approach)
    for( unsigned iter = 0; iter < _nIter; iter++ ) {
        if( rkStages == 2 ) solRK0 = _sol;

        for( unsigned rkStep = 0; rkStep < rkStages; ++rkStep ) {
            // --- Prepare Boundaries and temp variables
            IterPreprocessing( iter + rkStep );

            // --- Compute Fluxes ---
            FluxUpdate();

            // --- Finite Volume Update ---
            FVMUpdate( iter + rkStep );

            // --- Update Solution within Runge Kutta Stages
            _sol = _solNew;
        }
        // --- Iter Postprocessing ---
        IterPostprocessing( iter );

        // --- Wall time measurement
        duration  = std::chrono::high_resolution_clock::now() - start;
        _currTime = std::chrono::duration_cast<std::chrono::duration<double>>( duration ).count();

        // --- Runge Kutta Timestep ---
        if( rkStages == 2 ) RKUpdate( solRK0, _sol );

        // --- Write Output ---
        WriteVolumeOutput( iter );
        WriteScalarOutput( iter );

        // --- Update Scalar Fluxes
        _scalarFlux = _scalarFluxNew;

        // --- Print Output ---
        PrintScreenOutput( iter );
        PrintHistoryOutput( iter );
        PrintVolumeOutput( iter );
    }

    // --- Postprocessing ---

    DrawPostSolverOutput();
}

void SNSolverHPC::RKUpdate( VectorVector sol_0, VectorVector sol_rk ) {
#pragma omp parallel for
    for( unsigned i = 0; i < _nCells*_nSystem; ++i ) {
        _sol[i] = 0.5 * ( sol_0[i] + sol_rk[i] );
    }
}

void SNSolverHPC::PrintVolumeOutput() const { ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh ); }

void SNSolverHPC::PrintVolumeOutput( int currEnergy ) const {
    if( _settings->GetVolumeOutputFrequency() != 0 && currEnergy % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) {
        ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), _outputFields, _outputFieldNames, _mesh );
    }
    if( currEnergy == (int)_nIter - 1 ) {    // Last iteration write without suffix.
        ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh );
    }
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
            case MASS: _screenOutputFields[idx_field] = _problem->GetMass(); break;
            case ITER: _screenOutputFields[idx_field] = idx_iter; break;
            case WALL_TIME: _screenOutputFields[idx_field] = _currTime; break;
            case RMS_FLUX: _screenOutputFields[idx_field] = _problem->GetChangeRateFlux(); break;
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
            case CUR_OUTFLOW: _screenOutputFields[idx_field] = _problem->GetCurrentOutflow(); break;
            case TOTAL_OUTFLOW: _screenOutputFields[idx_field] = _problem->GetTotalOutflow(); break;
            case MAX_OUTFLOW: _screenOutputFields[idx_field] = _problem->GetMaxOrdinatewiseOutflow(); break;
            case CUR_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _problem->GetCurAbsorptionLattice(); break;
            case TOTAL_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _problem->GetTotalAbsorptionLattice(); break;
            case MAX_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _problem->GetMaxAbsorptionLattice(); break;
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
            case MASS: _historyOutputFields[idx_field] = _problem->GetMass(); break;
            case ITER: _historyOutputFields[idx_field] = idx_iter; break;
            case WALL_TIME: _historyOutputFields[idx_field] = _currTime; break;
            case RMS_FLUX: _historyOutputFields[idx_field] = _problem->GetChangeRateFlux(); break;
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
            case CUR_OUTFLOW: _historyOutputFields[idx_field] = _problem->GetCurrentOutflow(); break;
            case TOTAL_OUTFLOW: _historyOutputFields[idx_field] = _problem->GetTotalOutflow(); break;
            case MAX_OUTFLOW: _historyOutputFields[idx_field] = _problem->GetMaxOrdinatewiseOutflow(); break;
            case CUR_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _problem->GetCurAbsorptionLattice(); break;
            case TOTAL_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _problem->GetTotalAbsorptionLattice(); break;
            case MAX_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _problem->GetMaxAbsorptionLattice(); break;
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
    else if( idx_iter == _nEnergies - 1 ) {    // Always print last iteration
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
    log->info( "| The simulation will run for {} iterations.", _nEnergies );
    log->info( "| The spatial grid contains {} cells.", _nCells );
    if( _settings->GetSolverName() != PN_SOLVER && _settings->GetSolverName() != CSD_PN_SOLVER ) {
        log->info( "| The velocity grid contains {} points.", _quadrature->GetNq() );
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

void SNSolverHPC::SolverPreprocessing() {}

void SNSolverHPC::IterPostprocessing( unsigned /*idx_iter*/ ) {
    // --- Compute Quantities of interest for Volume and Screen Output ---
    ComputeScalarFlux();    // Needs to be called first is a solver function

    _problem->ComputeMass( _scalarFluxNew );
    _problem->ComputeChangeRateFlux( _scalarFlux, _scalarFluxNew );

    _problem->ComputeCurrentOutflow( _sol );
    _problem->ComputeTotalOutflow( _dT );
    _problem->ComputeMaxOrdinatewiseOutflow( _sol );

    if( _settings->GetProblemName() == PROBLEM_Lattice || _settings->GetProblemName() == PROBLEM_HalfLattice ) {
        _problem->ComputeCurrentAbsorptionLattice( _scalarFlux );
        _problem->ComputeTotalAbsorptionLattice( _dT );
        _problem->ComputeMaxAbsorptionLattice( _scalarFlux );
    }
    if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum || _settings->GetProblemName() == PROBLEM_QuarterHohlraum ) {
        _problem->ComputeCurrentAbsorptionHohlraum( _scalarFlux );    // Unify
        _problem->ComputeTotalAbsorptionHohlraum( _dT );              // Unify and parallelize
        _problem->ComputeCurrentProbeMoment( _sol );
        _problem->ComputeVarAbsorptionGreen( _scalarFlux );
        _problem->ComputeQOIsGreenProbingLine( _scalarFlux );
    }
}

