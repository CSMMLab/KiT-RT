#include "solvers/solverbase.h"
#include "common/config.h"
#include "common/globalconstants.h"
#include "common/io.h"
#include "common/mesh.h"
#include "fluxes/numericalflux.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"
#include "solvers/csdpnsolver.h"
#include "solvers/csdsnsolver.h"
#include "solvers/csdsolvertrafofp.h"
#include "solvers/csdsolvertrafofp2d.h"
#include "solvers/csdsolvertrafofpsh2d.h"
#include "toolboxes/textprocessingtoolbox.h"

#include "solvers/mnsolver.h"
#include "solvers/pnsolver.h"
#include "solvers/snsolver.h"

#include <mpi.h>

SolverBase::SolverBase( Config* settings ) {
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

    _interfaceMidPoints = _mesh->GetInterfaceMidPoints();

    _cellMidPoints = _mesh->GetCellMidPoints();

    // set time step or energy step
    _dE = ComputeTimeStep( _settings->GetCFL() );

    if( _settings->GetIsCSD() ) {
        // carefull: This gets overwritten by almost all subsolvers
        double minE = 5e-5;    // 2.231461e-01;    // 5e-5;
        double maxE = _settings->GetMaxEnergyCSD();
        _nEnergies  = std::ceil( ( maxE - minE ) / _dE );
        _energies   = blaze::linspace( _nEnergies, minE, maxE );
        //_energies = blaze::linspace( _nEnergies, maxE, minE );    // go backwards from biggest to smallest energy
    }
    else {    // Not CSD Solver
        _nEnergies = unsigned( settings->GetTEnd() / _dE );
        _energies  = blaze::linspace( _nEnergies, 0.0, settings->GetTEnd() );    // go upward from 0 to T_end
    }

    // setup problem  and store frequently used params
    _problem      = ProblemBase::Create( _settings, _mesh );
    _sol          = _problem->SetupIC();
    _solNew       = _sol;    // setup temporary sol variable
    unsigned nSys = _sol[0].size();
    _solDx        = VectorVector( _nCells, Vector( nSys, 0.0 ) );
    _solDy        = VectorVector( _nCells, Vector( nSys, 0.0 ) );
    _limiter      = VectorVector( _nCells, Vector( nSys, 0.0 ) );

    _sigmaT = _problem->GetTotalXS( _energies );
    _sigmaS = _problem->GetScatteringXS( _energies );
    _Q      = _problem->GetExternalSource( _energies );

    // setup numerical flux
    _g = NumericalFlux::Create();

    // boundary type
    _boundaryCells = _mesh->GetBoundaryTypes();

    // Solver Output
    _solverOutput.resize( _nCells );    // LEGACY! Only used for CSD SN

    PrepareScreenOutput();     // Screen Output
    PrepareHistoryOutput();    // History Output

    // initialize Helper Variables
    _fluxNew = Vector( _nCells, 0 );
    _flux    = Vector( _nCells, 0 );

    // write density
    _density = _problem->GetDensity( _mesh->GetCellMidPoints() );
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
        case CSD_SN_SOLVER: return new CSDSNSolver( settings );
        case CSD_SN_FOKKERPLANCK_TRAFO_SOLVER: return new CSDSolverTrafoFP( settings );
        case CSD_SN_FOKKERPLANCK_TRAFO_SOLVER_2D: return new CSDSolverTrafoFP2D( settings );
        case CSD_SN_FOKKERPLANCK_TRAFO_SH_SOLVER_2D: return new CSDSolverTrafoFPSH2D( settings );
        case CSD_PN_SOLVER: return new CSDPNSolver( settings );
        default: ErrorMessages::Error( "Creator for the chosen solver does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    }
    ErrorMessages::Error( "Creator for the chosen solver does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    return nullptr;    // This code is never reached. Just to disable compiler warnings.
}

void SolverBase::Solve() {
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

void SolverBase::PrintVolumeOutput() const { ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh ); }

void SolverBase::PrintVolumeOutput( int currEnergy ) const {
    if( _settings->GetVolumeOutputFrequency() != 0 && currEnergy % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) {
        ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), _outputFields, _outputFieldNames, _mesh );
    }
    if( currEnergy == (int)_maxIter - 1 ) {    // Last iteration write without suffix.
        ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh );
    }
}

// --- Helper ---
double SolverBase::ComputeTimeStep( double cfl ) const {
    double maxEdge = -1.0;
    for( unsigned j = 0; j < _nCells; j++ ) {
        for( unsigned l = 0; l < _normals[j].size(); l++ ) {
            double currentEdge = _areas[j] / norm( _normals[j][l] );
            if( currentEdge > maxEdge ) maxEdge = currentEdge;
        }
    }
    return cfl * maxEdge;
}

// --- High order reconstruction ---
void SolverBase::ComputeGradients( unsigned nSys ) {
    /*
    for( unsigned idx_cell = 0; idx_cell < _numCells; ++idx_cell ) {
        for( unsigned idx_sys = 0; idx_sys < nSys; ++idx_sys ) {
            _solDx[idx_cell][idx_sys] = 0.0;
            _solDy[idx_cell][idx_sys] = 0.0;
        }
        if( _boundaryCells[idx_cell] != BOUNDARY_TYPE::NONE ) continue;    // skip boundary cells

        for( unsigned idx_sys = 0; idx_sys < nSys; ++idx_sys ) {
            // compute derivative by summing over cell boundary using green gauss theorem and midpoint rule
            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); ++idx_nbr ) {
                _solDx[idx_cell][idx_sys] +=
                    0.5 * ( psi[idx_cell][idx_sys] + _sol[_cellNeighbors[idx_cell][idx_nbr]][idx_sys] ) * _normals[idx_cell][idx_nbr][0];
                _solDy[idx_cell][idx_sys] +=
                    0.5 * ( psi[idx_cell][idx_sys] + _sol[_cellNeighbors[idx_cell][idx_nbr]][idx_sys] ) * _normals[idx_cell][idx_nbr][1];
            }
            psiDerX[idx_cell][idx_sys] /= _cellAreas[idx_cell];
            psiDerY[idx_cell][idx_sys] /= _cellAreas[idx_cell];
        }
    }
    */
}

void SolverBase::ComputeLimiter() {
    /*
    double r    = 0.0;
    double eps  = 1e-10;
    double sign = 0.0;
    for( unsigned idx_cell = 0; idx_cell < _numCells; idx_cell++ ) {
        for( unsigned idx_sys = 0; idx_sys < nSys; idx_sys++ ) {
            if( _cellBoundaryTypes[idx_cell] != 2 ) {
                _limiter[idx_cell][idx_sys] = 0.0;    // turn to first order on boundaries
                continue;                             // skip computation
            }
            double minSol = _sol[idx_cell][idx_sys];
            double maxSol = _sol[idx_cell][idx_sys];
            Vector localLimiter( _numNodesPerCell, 0.0 );
            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); idx_nbr++ ) {
                // Compute ptswise max and minimum solultion values of current and neighbor cells
                unsigned glob_nbr = _cellNeighbors[idx_cell][idx_nbr];
                if( _sol[glob_nbr][idx_sys] > maxSol ) {
                    maxSol = _sol[glob_nbr][idx_sys];
                }
                if( sol[glob_nbr][idx_sys] < minSol ) {
                    minSol = _sol[glob_nbr][idx_sys];
                }
            }
            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); idx_nbr++ ) {
                // Compute value at interface midpoint, called gaussPt
                if( idx_cell == 1049 && idx_nbr == 1 ) {
                    // std::cout << "here\n";
                }
                double gaussPt = 0.0;
                double dy      = _solDy[idx_cell][idx_sys];
                double dx      = _solDx[idx_cell][idx_sys];
                double rijx    = _interfaceMidPoints[idx_cell][idx_nbr][0];
                double rijy    = _interfaceMidPoints[idx_cell][idx_nbr][1];
                double cmx     = _cellMidPoints[idx_cell][0];
                double cmy     = _cellMidPoints[idx_cell][1];
                double curSol  = _sol[idx_cell][idx_sys];

                gaussPt = _solDx[idx_cell][idx_sys] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[idx_cell][0] ) +
                          _solDy[idx_cell][idx_sys] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[idx_cell][1] );
                // Compute limiter input
                if( gaussPt < 0.0 )
                    sign = -1.0;
                else
                    sign = 1.0;

                if( gaussPt > 0.0 ) {
                    double t1 = maxSol - _sol[idx_cell][idx_sys];
                    r         = ( maxSol - _sol[idx_cell][idx_sys] ) / ( sign * ( std::abs( gaussPt ) + eps ) );
                }
                else if( gaussPt < 0.0 ) {
                    double t1 = minSol - _sol[idx_cell][idx_sys];
                    r         = ( minSol - _sol[idx_cell][idx_sys] ) / ( sign * ( std::abs( gaussPt ) + eps ) );
                }
                else
                    r = 1.0;
                if( r < 0 ) {
                    std::cout << "here\n";
                }
                localLimiter[idx_nbr] = std::min( r, 1.0 );    // LimiterBarthJespersen( r );
            }
            // get smallest limiter
            limiter[idx_cell][idx_sys] = localLimiter[0];
            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); idx_nbr++ ) {
                if( localLimiter[idx_nbr] < limiter[idx_cell][idx_sys] ) limiter[idx_cell][idx_sys] = localLimiter[idx_nbr];
            }
            // check maximum principle
            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); idx_nbr++ ) {
                double currLim = limiter[idx_cell][idx_sys];
                double dy      = _solDy[idx_cell][idx_sys];
                double dx      = _solDx[idx_cell][idx_sys];
                double rijx    = _interfaceMidPoints[idx_cell][idx_nbr][0];
                double rijy    = _interfaceMidPoints[idx_cell][idx_nbr][1];
                double cmx     = _cellMidPoints[idx_cell][0];
                double cmy     = _cellMidPoints[idx_cell][1];
                double curSol  = _sol[idx_cell][idx_sys];
                double gaussPt = _solDx[idx_cell][idx_sys] * ( _interfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[idx_cell][0] ) +
                                 _solDy[idx_cell][idx_sys] * ( _interfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[idx_cell][1] );

                double psiL  = _sol[idx_cell][idx_sys] + currLim * gaussPt;
                double psiL2 = curSol + currLim * ( dx * ( rijx - cmx ) + dy * ( rijy - cmy ) );

                if( psiL > maxSol ) {
                    std::cout << "here\n";
                }
                if( psiL < minSol ) {
                    std::cout << "here\n";
                }
            }
        }
    }
    */
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

            default: ErrorMessages::Error( "Screen output field not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void SolverBase::WriteScalarOutput( unsigned iteration ) {

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

            case RMS_FLUX: _screenOutputFields[idx_field] = blaze::l2Norm( _fluxNew - _flux ); break;

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
    _flux = _fluxNew;
}

void SolverBase::PrintScreenOutput( unsigned iteration ) {
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

            default: ErrorMessages::Error( "History output field not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void SolverBase::PrintHistoryOutput( unsigned iteration ) {
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

void SolverBase::DrawPreSolverOutput() {
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
        log->info( "| Postprocessing screen output goes here." );
        log->info( "--------------------------- Solver Finished ----------------------------" );
    }
}

void SolverBase::SolverPreprocessing() {}
