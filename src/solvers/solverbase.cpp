#include "solvers/solverbase.hpp"
#include "common/config.hpp"
#include "common/globalconstants.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "fluxes/numericalflux.hpp"
#include "problems/problembase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "solvers/csdmnsolver.hpp"
#include "solvers/csdpnsolver.hpp"
#include "solvers/csdsnsolver.hpp"
#include "solvers/mnsolver.hpp"
#include "solvers/mnsolver_normalized.hpp"
#include "solvers/pnsolver.hpp"
#include "solvers/snsolver.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include <mpi.h>

SolverBase::SolverBase( Config* settings ) {
    _settings = settings;
    // build mesh and store  and store frequently used params
    _mesh = LoadSU2MeshFromFile( settings );

    _areas     = _mesh->GetCellAreas();
    _neighbors = _mesh->GetNeighbours();
    _normals   = _mesh->GetNormals();
    _nCells    = _mesh->GetNumCells();
    _settings->SetNCells( _nCells );

    // build quadrature object and store frequently used params
    _quadrature = QuadratureBase::Create( settings );
    _nq         = _quadrature->GetNq();
    _settings->SetNQuadPoints( _nq );

    // build slope related params
    _reconstructor      = new Reconstructor( settings );    // Not used!
    _reconsOrder        = _reconstructor->GetReconsOrder();
    _interfaceMidPoints = _mesh->GetInterfaceMidPoints();

    _cellMidPoints = _mesh->GetCellMidPoints();

    // set time step or energy step
    _dT = ComputeTimeStep( _settings->GetCFL() );

    if( _settings->GetIsCSD() ) {
        // carefull: This gets overwritten by almost all subsolvers
        double minE = 5e-5;    // 2.231461e-01;    // 5e-5;
        double maxE = _settings->GetMaxEnergyCSD();
        _nEnergies  = std::ceil( ( maxE - minE ) / _dT );
        _energies   = blaze::linspace( _nEnergies, minE, maxE );
    }
    else {                                                     // Not CSD Solver
        _nEnergies = unsigned( settings->GetTEnd() / _dT );    // redundancy with nIter (<-Better) ?
        _energies  = 0;    // blaze::linspace( _nEnergies, 0.0, settings->GetTEnd() );    // go upward from 0 to T_end =>Not needed
    }

    // Adjust maxIter, depending if we have a normal run or a csd Run
    _nIter = _nEnergies;
    if( _settings->GetIsCSD() ) {
        _nIter = _nEnergies - 1;    // Since CSD does not go the last energy step
    }

    // setup problem  and store frequently used params

    _problem = ProblemBase::Create( _settings, _mesh, _quadrature );
    _sol     = _problem->SetupIC();

    _solNew = _sol;    // setup temporary sol variable
    if( !_settings->GetIsCSD() ) {
        _sigmaT = _problem->GetTotalXS( _energies );
        _sigmaS = _problem->GetScatteringXS( _energies );
        _Q      = _problem->GetExternalSource( _energies );
    }

    // setup numerical flux
    _g = NumericalFluxBase::Create();

    // boundary type
    _boundaryCells = _mesh->GetBoundaryTypes();

    // Solver Output
    _solverOutput.resize( _nCells );    // LEGACY! Only used for CSD SN

    PrepareScreenOutput();     // Screen Output
    PrepareHistoryOutput();    // History Output

    // initialize Helper Variables
    _scalarFluxNew = Vector( _nCells, 0 );
    _scalarFlux    = Vector( _nCells, 0 );

    // write density
    _density = _problem->GetDensity( _mesh->GetCellMidPoints() );

    // initialize QOI helper variables
    _curMaxOrdinateOutflow = 0.0;
    _curScalarOutflow      = 0.0;
    _totalScalarOutflow    = 0.0;
    _mass                  = 0.0;
    _changeRateFlux        = 0.0;

    // Hardcoded for the symmetric hohlraum testcase (experimental, needs refactoring)
    if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) {
        _probingCells = {
            _mesh->GetCellOfKoordinate( -0.4, 0. ),
            _mesh->GetCellOfKoordinate( 0.4, 0. ),
            _mesh->GetCellOfKoordinate( 0., -0.6 ),
            _mesh->GetCellOfKoordinate( 0., 0.6 ),
        };
        _probingMoments = VectorVector( 4, Vector( 3, 0.0 ) );
    }
    else {
        _probingCells   = { 0 };    // Dummy needs refactoring
        _probingMoments = VectorVector( 4, Vector( 3, 0.0 ) );
    }
}

SolverBase::~SolverBase() {
    delete _quadrature;
    delete _mesh;
    delete _problem;
    delete _reconstructor;
    delete _g;
}

SolverBase* SolverBase::Create( Config* settings ) {
    switch( settings->GetSolverName() ) {
        case SN_SOLVER: return new SNSolver( settings );
        case PN_SOLVER: return new PNSolver( settings );
        case MN_SOLVER: return new MNSolver( settings );
        case MN_SOLVER_NORMALIZED: return new MNSolverNormalized( settings );
        case CSD_SN_SOLVER: return new CSDSNSolver( settings );
        case CSD_PN_SOLVER: return new CSDPNSolver( settings );
        case CSD_MN_SOLVER: return new CSDMNSolver( settings );

        default: ErrorMessages::Error( "Creator for the chosen solver does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    }
    ErrorMessages::Error( "Creator for the chosen solver does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    return nullptr;    // This code is never reached. Just to disable compiler warnings.
}

void SolverBase::Solve() {

    // --- Preprocessing ---

    PrepareVolumeOutput();

    DrawPreSolverOutput();

    // Preprocessing before first pseudo time step
    SolverPreprocessing();
    unsigned rkStages = _settings->GetRKStages();
    // Create Backup solution for Runge Kutta
    VectorVector solRK0 = _sol;

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

void SolverBase::RKUpdate( VectorVector sol_0, VectorVector sol_rk ) {
#pragma omp parallel for
    for( unsigned i = 0; i < _nCells; ++i ) {
        _sol[i] = 0.5 * ( sol_0[i] + sol_rk[i] );
    }
}

void SolverBase::PrintVolumeOutput() const { ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh ); }

void SolverBase::PrintVolumeOutput( int currEnergy ) const {
    if( _settings->GetVolumeOutputFrequency() != 0 && currEnergy % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) {
        ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), _outputFields, _outputFieldNames, _mesh );
    }
    if( currEnergy == (int)_nIter - 1 ) {    // Last iteration write without suffix.
        ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh );
    }
}

// --- Helper ---
double SolverBase::ComputeTimeStep( double cfl ) const {
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
    double maxEdge = -1.0;
    for( unsigned j = 0; j < _nCells; j++ ) {
        for( unsigned l = 0; l < _normals[j].size(); l++ ) {
            double currentEdge = _areas[j] / norm( _normals[j][l] );
            if( currentEdge > maxEdge ) maxEdge = currentEdge;
        }
    }
    return cfl * maxEdge;
}

// --- IO ----
void SolverBase::PrepareScreenOutput() {
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
                idx_field++;
                _screenOutputFieldNames[idx_field] = "Probe 3 u_0";
                idx_field++;
                _screenOutputFieldNames[idx_field] = "Probe 4 u_0";
                break;
            case VAR_ABSORPTION_GREEN: _screenOutputFieldNames[idx_field] = "Var. absorption green"; break;
            default: ErrorMessages::Error( "Screen output field not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void SolverBase::WriteScalarOutput( unsigned idx_iter ) {

    unsigned nFields = (unsigned)_settings->GetNScreenOutput();
    double mass      = 0.0;

    // -- Screen Output
    for( unsigned idx_field = 0; idx_field < nFields; idx_field++ ) {
        // Prepare all Output Fields per group
        // Different procedure, depending on the Group...
        switch( _settings->GetScreenOutput()[idx_field] ) {
            case MASS: _screenOutputFields[idx_field] = _mass; break;
            case ITER: _screenOutputFields[idx_field] = idx_iter; break;
            case RMS_FLUX: _screenOutputFields[idx_field] = _changeRateFlux; break;
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
            case CUR_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _problem->GetCurAbsorptionLattice(); break;
            case TOTAL_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _problem->GetTotalAbsorptionLattice(); break;
            case MAX_PARTICLE_ABSORPTION: _screenOutputFields[idx_field] = _problem->GetMaxAbsorptionLattice(); break;
            case TOTAL_PARTICLE_ABSORPTION_CENTER: _screenOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumCenter(); break;
            case TOTAL_PARTICLE_ABSORPTION_VERTICAL: _screenOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumVertical(); break;
            case TOTAL_PARTICLE_ABSORPTION_HORIZONTAL: _screenOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumHorizontal(); break;
            case PROBE_MOMENT_TIME_TRACE:
                for( unsigned i = 0; i < 4; i++ ) {
                    _screenOutputFields[idx_field] = _probingMoments[i][0];
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
            case RMS_FLUX: _historyOutputFields[idx_field] = _changeRateFlux; break;
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
            case CUR_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _problem->GetCurAbsorptionLattice(); break;
            case TOTAL_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _problem->GetTotalAbsorptionLattice(); break;
            case MAX_PARTICLE_ABSORPTION: _historyOutputFields[idx_field] = _problem->GetMaxAbsorptionLattice(); break;
            case TOTAL_PARTICLE_ABSORPTION_CENTER: _historyOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumCenter(); break;
            case TOTAL_PARTICLE_ABSORPTION_VERTICAL: _historyOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumVertical(); break;
            case TOTAL_PARTICLE_ABSORPTION_HORIZONTAL: _historyOutputFields[idx_field] = _problem->GetTotalAbsorptionHohlraumHorizontal(); break;
            case PROBE_MOMENT_TIME_TRACE:
                for( unsigned i = 0; i < 4; i++ ) {
                    for( unsigned j = 0; j < 3; j++ ) {
                        _historyOutputFields[idx_field] = _probingMoments[i][j];
                        idx_field++;
                    }
                }
                idx_field--;
                break;
            case VAR_ABSORPTION_GREEN: _historyOutputFields[idx_field] = _problem->GetVarAbsorptionHohlraumGreen(); break;

            default: ErrorMessages::Error( "History output group not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void SolverBase::PrintScreenOutput( unsigned idx_iter ) {
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
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
                                                        CUR_OUTFLOW,
                                                        TOTAL_OUTFLOW,
                                                        MAX_OUTFLOW,
                                                        CUR_PARTICLE_ABSORPTION,
                                                        TOTAL_PARTICLE_ABSORPTION,
                                                        MAX_PARTICLE_ABSORPTION,
                                                        TOTAL_PARTICLE_ABSORPTION_CENTER,
                                                        TOTAL_PARTICLE_ABSORPTION_VERTICAL,
                                                        TOTAL_PARTICLE_ABSORPTION_HORIZONTAL,
                                                        PROBE_MOMENT_TIME_TRACE };
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
            ss << _screenOutputFields[idx_field];
            tmp = ss.str();
            tmp.erase( std::remove( tmp.begin(), tmp.end(), '+' ), tmp.end() );    // removing the '+' sign
        }

        if( strLen > tmp.size() )    // Padding
            tmp.insert( 0, strLen - tmp.size(), paddingChar );
        else if( strLen < tmp.size() )    // Cutting
            tmp.resize( strLen );

        lineToPrint += tmp + " |";
    }
    if( rank == 0 ) {
        if( _settings->GetScreenOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetScreenOutputFrequency() == 0 ) {
            log->info( lineToPrint );
        }
        else if( idx_iter == _nIter - 1 ) {    // Always print last iteration
            log->info( lineToPrint );
        }
    }
}

void SolverBase::PrepareHistoryOutput() {
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
            case RMS_FLUX: _historyOutputFieldNames[idx_field] = "RMS_flux"; break;
            case VTK_OUTPUT: _historyOutputFieldNames[idx_field] = "VTK_out"; break;
            case CSV_OUTPUT: _historyOutputFieldNames[idx_field] = "CSV_out"; break;
            case CUR_OUTFLOW: _historyOutputFieldNames[idx_field] = "Final_time_outflow"; break;
            case TOTAL_OUTFLOW: _historyOutputFieldNames[idx_field] = "Cumulated_outflow"; break;
            case MAX_OUTFLOW: _historyOutputFieldNames[idx_field] = "Max_outflow"; break;
            case CUR_PARTICLE_ABSORPTION: _historyOutputFieldNames[idx_field] = "Final_time_absorption"; break;
            case TOTAL_PARTICLE_ABSORPTION: _historyOutputFieldNames[idx_field] = "Cumulated_absorption"; break;
            case MAX_PARTICLE_ABSORPTION: _historyOutputFieldNames[idx_field] = "Max_absorption"; break;
            case TOTAL_PARTICLE_ABSORPTION_CENTER: _historyOutputFieldNames[idx_field] = "Cumulated_absorption_center"; break;
            case TOTAL_PARTICLE_ABSORPTION_VERTICAL: _historyOutputFieldNames[idx_field] = "Cumulated_absorption_vertical_wall"; break;
            case TOTAL_PARTICLE_ABSORPTION_HORIZONTAL: _historyOutputFieldNames[idx_field] = "Cumulated_absorption_horizontal_wall"; break;
            case PROBE_MOMENT_TIME_TRACE:
                for( unsigned i = 0; i < 4; i++ ) {
                    for( unsigned j = 0; j < 3; j++ ) {
                        _historyOutputFieldNames[idx_field] = "Probe " + std::to_string( i ) + " u_" + std::to_string( j );
                        idx_field++;
                    }
                }
                idx_field--;
                break;
            case VAR_ABSORPTION_GREEN: _historyOutputFieldNames[idx_field] = "Var. absorption green"; break;
            default: ErrorMessages::Error( "History output field not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void SolverBase::PrintHistoryOutput( unsigned idx_iter ) {
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    auto log = spdlog::get( "tabular" );

    // assemble the line to print
    std::string lineToPrint = "";
    std::string tmp;
    for( int idx_field = 0; idx_field < _settings->GetNScreenOutput() - 1; idx_field++ ) {
        tmp = std::to_string( _screenOutputFields[idx_field] );
        lineToPrint += tmp + ",";
    }
    tmp = std::to_string( _screenOutputFields[_settings->GetNScreenOutput() - 1] );
    lineToPrint += tmp;    // Last element without comma

    if( rank == 0 ) {
        if( _settings->GetHistoryOutputFrequency() != 0 && idx_iter % (unsigned)_settings->GetHistoryOutputFrequency() == 0 ) {
            log->info( lineToPrint );
        }
        else if( idx_iter == _nEnergies - 1 ) {    // Always print last iteration
            log->info( lineToPrint );
        }
    }
}

void SolverBase::DrawPreSolverOutput() {
    // MPI
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // Logger
    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );

    std::string hLine = "--";

    if( rank == 0 ) {
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
}

void SolverBase::DrawPostSolverOutput() {
    // MPI
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // Logger
    auto log = spdlog::get( "event" );

    std::string hLine = "--";

    if( rank == 0 ) {
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
}

void SolverBase::SolverPreprocessing() {}

void SolverBase::GetTotalOutflow() { _totalScalarOutflow += _curScalarOutflow * _dT; }

void SolverBase::GetMass() {
    _mass = 0.0;
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        _mass += _scalarFluxNew[idx_cell] * _areas[idx_cell];
    }
}

void SolverBase::GetChangeRateFlux() { _changeRateFlux = blaze::l2Norm( _scalarFluxNew - _scalarFlux ); }

void SolverBase::IterPostprocessing( unsigned idx_iter ) {
    // --- Compute Quantities of interest for Volume and Screen Output ---
    ComputeScalarFlux();    // Needs to be called first

    GetMass();
    GetCurrentOutflow();
    GetTotalOutflow();

    GetMaxOrdinatewiseOutflow();

    if( _settings->GetProblemName() == PROBLEM_Lattice ) {
        _problem->ComputeCurrentAbsorptionLattice( _scalarFlux );
        _problem->ComputeTotalAbsorptionLattice( _dT );
        _problem->ComputeMaxAbsorptionLattice( _scalarFlux );
    }
    if( _settings->GetProblemName() == PROBLEM_SymmetricHohlraum ) {
        _problem->ComputeCurrentAbsorptionHohlraum( _scalarFlux );
        _problem->ComputeTotalAbsorptionHohlraum( _dT );
        ComputeCurrentProbeMoment();
        _problem->ComputeVarAbsorptionGreen( _scalarFlux );
    }
}
