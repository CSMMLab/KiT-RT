#include "solvers/solverbase.h"
#include "common/config.h"
#include "common/globalconstants.h"
#include "common/io.h"
#include "common/mesh.h"
#include "fluxes/numericalflux.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"
#include "solvers/csdsnsolver.h"
#include "solvers/csdsolvertrafofp.h"
#include "solvers/csdsolvertrafofp2d.h"
#include "solvers/csdsolvertrafofpsh2d.h"
#include "toolboxes/textprocessingtoolbox.h"

#include "solvers/mnsolver.h"
#include "solvers/pnsolver.h"
#include "solvers/snsolver.h"

#include <mpi.h>

Solver::Solver( Config* settings ) {
    _settings = settings;

    // @TODO save parameters from settings class

    // build mesh and store  and store frequently used params
    _mesh      = LoadSU2MeshFromFile( settings );
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
    _reconstructor = new Reconstructor( settings );
    _reconsOrder   = _reconstructor->GetReconsOrder();

    auto nodes = _mesh->GetNodes();
    auto cells = _mesh->GetCells();
    std::vector<std::vector<Vector>> interfaceMidPoints( _nCells, std::vector<Vector>( _mesh->GetNumNodesPerCell(), Vector( 2, 1e-10 ) ) );
    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
        for( unsigned k = 0; k < _mesh->GetDim(); ++k ) {
            for( unsigned j = 0; j < _neighbors[idx_cell].size() - 1; ++j ) {
                interfaceMidPoints[idx_cell][j][k] = 0.5 * ( nodes[cells[idx_cell][j]][k] + nodes[cells[idx_cell][j + 1]][k] );
            }
            interfaceMidPoints[idx_cell][_neighbors[idx_cell].size() - 1][k] =
                0.5 * ( nodes[cells[idx_cell][_neighbors[idx_cell].size() - 1]][k] + nodes[cells[idx_cell][0]][k] );
        }
    }
    _interfaceMidPoints = interfaceMidPoints;
    _cellMidPoints      = _mesh->GetCellMidPoints();

    _psiDx = VectorVector( _nCells, Vector( _nq, 0.0 ) );
    _psiDy = VectorVector( _nCells, Vector( _nq, 0.0 ) );

    // set time step
    _dE        = ComputeTimeStep( settings->GetCFL() );
    _nEnergies = unsigned( settings->GetTEnd() / _dE );
    _energies.resize( _nEnergies );
    for( unsigned i = 0; i < _nEnergies; ++i ) _energies[i] = ( i + 1 ) * _dE;

    // setup problem  and store frequently used params
    _problem      = ProblemBase::Create( _settings, _mesh );
    _sol          = _problem->SetupIC();
    auto cellMids = _mesh->GetCellMidPoints();
    // double enterPositionX = 0.5;    // 0.0;
    // double enterPositionY = 0.5;
    // auto boundaryCells    = _mesh->GetBoundaryTypes();
    // Case 1: Ingoing radiation in just one cell
    // find cell that best matches enter position
    // double dist    = 1000.0;
    _boundaryCells = _mesh->GetBoundaryTypes();
    // unsigned indexSource = 0;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        // if( boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];
        if( x >= 0.49 && x <= 0.5 && y >= 0.49 && y <= 0.5 ) {
            //_boundaryCells[j] = DIRICHLET;
        }
        //}
    }
    _solNew = _sol;    // setup temporary sol variable

    _sigmaT = _problem->GetTotalXS( _energies );
    _sigmaS = _problem->GetScatteringXS( _energies );
    _Q      = _problem->GetExternalSource( _energies );

    // setup numerical flux
    _g = NumericalFlux::Create();

    // boundary type
    //_boundaryCells              = _mesh->GetBoundaryTypes();
    //_boundaryCells[indexSource] = DIRICHLET;

    // Solver Output
    _solverOutput.resize( _nCells );    // LEGACY! Only used for CSD SN

    PrepareScreenOutput();     // Screen Output
    PrepareHistoryOutput();    // History Output

    // initialize Helper Variables
    _fluxNew = Vector( _nCells, 0 );
    _flux    = Vector( _nCells, 0 );

    // write density
    _density = _problem->GetDensity( _mesh->GetCellMidPoints() );
    //_density = std::vector( _mesh->GetCellMidPoints().size(), 0.0 );
}

Solver::~Solver() {
    delete _quadrature;
    delete _mesh;
    delete _problem;
}

Solver* Solver::Create( Config* settings ) {
    switch( settings->GetSolverName() ) {
        case SN_SOLVER: return new SNSolver( settings );

        case PN_SOLVER: return new PNSolver( settings );
        case MN_SOLVER: return new MNSolver( settings );
        case CSD_SN_SOLVER: return new CSDSNSolver( settings );
        case CSD_SN_FOKKERPLANCK_TRAFO_SOLVER: return new CSDSolverTrafoFP( settings );
        case CSD_SN_FOKKERPLANCK_TRAFO_SOLVER_2D: return new CSDSolverTrafoFP2D( settings );
        case CSD_SN_FOKKERPLANCK_TRAFO_SH_SOLVER_2D: return new CSDSolverTrafoFPSH2D( settings );
        default: ErrorMessages::Error( "Creator for the chosen solver does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    }
    ErrorMessages::Error( "Creator for the chosen solver does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    return nullptr;    // This code is never reached. Just to disable compiler warnings.
}

void Solver::Solve() {

    // --- Preprocessing ---

    PrepareVolumeOutput();

    DrawPreSolverOutput();

    // Adjust maxIter, depending if we have a normal run or a csd Run
    _maxIter = _nEnergies;
    if( _settings->GetIsCSD() ) {
        _maxIter = _nEnergies - 1;    // Since CSD does not go the last energy step
    }

    // Preprocessing before first pseudo time step
    SolverPreprocessing();

    // Loop over energies (pseudo-time of continuous slowing down approach)
    for( unsigned iter = 0; iter < _maxIter; iter++ ) {

        // --- Prepare Boundaries and temp variables
        IterPreprocessing( iter );

        // --- Compute Fluxes ---
        FluxUpdate();

        // --- Finite Volume Update ---
        FVMUpdate( iter );

        // --- Iter Postprocessing ---
        IterPostprocessing( iter );

        // --- Solver Output ---
        WriteVolumeOutput( iter );
        WriteScalarOutput( iter );
        PrintScreenOutput( iter );
        PrintHistoryOutput( iter );
        PrintVolumeOutput( iter );
    }

    // --- Postprocessing ---

    DrawPostSolverOutput();
}

void Solver::PrintVolumeOutput() const { ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh ); }

void Solver::PrintVolumeOutput( int currEnergy ) const {
    if( _settings->GetVolumeOutputFrequency() != 0 && currEnergy % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) {
        ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), _outputFields, _outputFieldNames, _mesh );
    }
    if( currEnergy == (int)_maxIter - 1 ) {    // Last iteration write without suffix.
        ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh );
    }
}

// --- Helper ---
double Solver::ComputeTimeStep( double cfl ) const {
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
void Solver::PrepareScreenOutput() {
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

            default: ErrorMessages::Error( "Screen output field not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void Solver::WriteScalarOutput( unsigned iteration ) {

    unsigned nFields = (unsigned)_settings->GetNScreenOutput();
    double mass      = 0.0;

    // -- Screen Output
    for( unsigned idx_field = 0; idx_field < nFields; idx_field++ ) {
        // Prepare all Output Fields per group
        // Different procedure, depending on the Group...
        switch( _settings->GetScreenOutput()[idx_field] ) {
            case MASS:
                for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                    mass += _fluxNew[idx_cell] * _areas[idx_cell];
                }
                _screenOutputFields[idx_field] = mass;
                break;

            case ITER: _screenOutputFields[idx_field] = iteration; break;

            case RMS_FLUX:
                _screenOutputFields[idx_field] = blaze::l2Norm( _fluxNew - _flux );
                _flux                          = _fluxNew;
                break;

            case VTK_OUTPUT:
                _screenOutputFields[idx_field] = 0;
                if( ( _settings->GetVolumeOutputFrequency() != 0 && iteration % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
                    ( iteration == _maxIter - 1 ) /* need sol at last iteration */ ) {
                    _screenOutputFields[idx_field] = 1;
                }
                break;

            case CSV_OUTPUT:
                _screenOutputFields[idx_field] = 0;
                if( ( _settings->GetHistoryOutputFrequency() != 0 && iteration % (unsigned)_settings->GetHistoryOutputFrequency() == 0 ) ||
                    ( iteration == _maxIter - 1 ) /* need sol at last iteration */ ) {
                    _screenOutputFields[idx_field] = 1;
                }
                break;

            default: ErrorMessages::Error( "Screen output group not defined!", CURRENT_FUNCTION ); break;
        }
    }

    // --- History output ---
    nFields = (unsigned)_settings->GetNHistoryOutput();

    std::vector<SCALAR_OUTPUT> screenOutputFields = _settings->GetScreenOutput();
    for( unsigned idx_field = 0; idx_field < nFields; idx_field++ ) {

        // Check first, if the field was already filled by screenoutput writer!
        std::vector<SCALAR_OUTPUT>::iterator itScreenOutput =
            std::find( screenOutputFields.begin(), screenOutputFields.end(), _settings->GetHistoryOutput()[idx_field] );

        // Prepare all Output Fields per group
        // Different procedure, depending on the Group...
        switch( _settings->GetHistoryOutput()[idx_field] ) {
            case MASS:
                if( screenOutputFields.end() == itScreenOutput ) {
                    for( unsigned idx_cell = 0; idx_cell < _nCells; ++idx_cell ) {
                        mass += _fluxNew[idx_cell] * _areas[idx_cell];
                    }
                    _historyOutputFields[idx_field] = mass;
                }
                else {
                    _historyOutputFields[idx_field] = *itScreenOutput;
                }
                break;

            case ITER: _historyOutputFields[idx_field] = iteration; break;

            case RMS_FLUX:
                if( screenOutputFields.end() == itScreenOutput ) {
                    _screenOutputFields[idx_field] = blaze::l2Norm( _fluxNew - _flux );
                    _flux                          = _fluxNew;
                }
                else {
                    _historyOutputFields[idx_field] = *itScreenOutput;
                }
                break;

            case VTK_OUTPUT:
                _historyOutputFields[idx_field] = 0;
                if( ( _settings->GetVolumeOutputFrequency() != 0 && iteration % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) ||
                    ( iteration == _maxIter - 1 ) /* need sol at last iteration */ ) {
                    _historyOutputFields[idx_field] = 1;
                }
                break;

            case CSV_OUTPUT:
                _historyOutputFields[idx_field] = 0;
                if( ( _settings->GetHistoryOutputFrequency() != 0 && iteration % (unsigned)_settings->GetHistoryOutputFrequency() == 0 ) ||
                    ( iteration == _maxIter - 1 ) /* need sol at last iteration */ ) {
                    _historyOutputFields[idx_field] = 1;
                }
                break;

            default: ErrorMessages::Error( "History output group not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void Solver::PrintScreenOutput( unsigned iteration ) {
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    auto log = spdlog::get( "event" );

    unsigned strLen  = 10;    // max width of one column
    char paddingChar = ' ';

    // assemble the line to print
    std::string lineToPrint = "| ";
    std::string tmp;
    for( unsigned idx_field = 0; idx_field < _settings->GetNScreenOutput(); idx_field++ ) {
        tmp = std::to_string( _screenOutputFields[idx_field] );

        // Format outputs correctly
        std::vector<SCALAR_OUTPUT> integerFields    = { ITER };
        std::vector<SCALAR_OUTPUT> scientificFields = { RMS_FLUX, MASS };
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
        if( _settings->GetScreenOutputFrequency() != 0 && iteration % (unsigned)_settings->GetScreenOutputFrequency() == 0 ) {
            log->info( lineToPrint );
        }
        else if( iteration == _maxIter - 1 ) {    // Always print last iteration
            log->info( lineToPrint );
        }
    }
}

void Solver::PrepareHistoryOutput() {
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

            default: ErrorMessages::Error( "History output field not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void Solver::PrintHistoryOutput( unsigned iteration ) {
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
        if( _settings->GetHistoryOutputFrequency() != 0 && iteration % (unsigned)_settings->GetHistoryOutputFrequency() == 0 ) {
            log->info( lineToPrint );
        }
        else if( iteration == _nEnergies - 1 ) {    // Always print last iteration
            log->info( lineToPrint );
        }
    }
}

void Solver::DrawPreSolverOutput() {
    // MPI
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // Logger
    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );

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
        log->info( "---------------------------- Solver Starts -----------------------------" );
        log->info( "| The simulation will run for {} iterations.", _nEnergies );
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

void Solver::DrawPostSolverOutput() {
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
        log->info( "| Postprocessing screen output goes here." );
        log->info( "--------------------------- Solver Finished ----------------------------" );
    }
}

void Solver::SolverPreprocessing() {}
