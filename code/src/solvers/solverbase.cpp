#include "solvers/solverbase.h"
#include "common/config.h"
#include "common/globalconstants.h"
#include "common/io.h"
#include "common/mesh.h"
#include "fluxes/numericalflux.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"
#include "solvers/csdsnsolver.h"
#include "solvers/mnsolver.h"
#include "solvers/pnsolver.h"
#include "solvers/snsolver.h"

Solver::Solver( Config* settings ) : _settings( settings ) {
    // @TODO save parameters from settings class

    // build mesh and store  and store frequently used params
    _mesh      = LoadSU2MeshFromFile( settings );
    _areas     = _mesh->GetCellAreas();
    _neighbors = _mesh->GetNeighbours();
    _normals   = _mesh->GetNormals();
    _nCells    = _mesh->GetNumCells();
    _settings->SetNCells( _nCells );

    // build quadrature object and store frequently used params
    _quadrature = QuadratureBase::CreateQuadrature( settings );
    _nq         = _quadrature->GetNq();
    _settings->SetNQuadPoints( _nq );

    // set time step
    _dE        = ComputeTimeStep( settings->GetCFL() );
    _nEnergies = unsigned( settings->GetTEnd() / _dE );
    _energies.resize( _nEnergies );
    for( unsigned i = 0; i < _nEnergies; ++i ) _energies[i] = ( i + 1 ) * _dE;

    // setup problem  and store frequently used params
    _problem = ProblemBase::Create( _settings, _mesh );
    _sol     = _problem->SetupIC();
    _solNew  = _sol;    // setup temporary sol variable

    //_s       = _problem->GetStoppingPower( _energies );
    _sigmaT = _problem->GetTotalXS( _energies );
    _sigmaS = _problem->GetScatteringXS( _energies );
    _Q      = _problem->GetExternalSource( _energies );

    // setup numerical flux
    _g = NumericalFlux::Create();

    // boundary type
    _boundaryCells = _mesh->GetBoundaryTypes();

    // Solver Output
    _solverOutput.resize( _nCells );    // LEGACY! Only used for CSD SN
    // Screen Output
    PrepareScreenOutputFields();

    // initialize Helper Variables
    _fluxNew = Vector( _nCells, 0 );
    _flux    = Vector( _nCells, 0 );
}

Solver::~Solver() {
    delete _quadrature;
    delete _mesh;
    delete _problem;
}

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

Solver* Solver::Create( Config* settings ) {
    switch( settings->GetSolverName() ) {
        case SN_SOLVER: return new SNSolver( settings );
        case CSD_SN_SOLVER: return new CSDSNSolver( settings );
        case PN_SOLVER: return new PNSolver( settings );
        case MN_SOLVER: return new MNSolver( settings );
        default: return new SNSolver( settings );
    }
}

void Solver::Save() const { ExportVTK( _settings->GetOutputFile(), _outputFields, _outputFieldNames, _mesh ); }

void Solver::Save( int currEnergy ) const {
    if( _settings->GetVolumeOutputFrequency() != 0 && currEnergy % (unsigned)_settings->GetVolumeOutputFrequency() == 0 ) {
        ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( currEnergy ), _outputFields, _outputFieldNames, _mesh );
    }
}

void Solver::PrepareScreenOutputFields() {
    unsigned nFields = (unsigned)_settings->GetNScreenOutput();

    _screenOutputFieldNames.resize( nFields );
    _screenOutputFields.resize( nFields );

    // Prepare all output Fields ==> Specified in option SCREEN_OUTPUT
    for( unsigned idx_field = 0; idx_field < nFields; idx_field++ ) {
        // Prepare all Output Fields per group

        // Different procedure, depending on the Group...
        switch( _settings->GetScreenOutput()[idx_field] ) {
            case MASS: _screenOutputFieldNames[idx_field] = "mass"; break;

            case ITER: _screenOutputFieldNames[idx_field] = "iter"; break;

            case RMS_FLUX: _screenOutputFieldNames[idx_field] = "RMS_flux"; break;

            default: ErrorMessages::Error( "Screen Output Group not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void Solver::WriteScreenOutputFields( unsigned idx_pseudoTime ) {

    unsigned nFields = (unsigned)_settings->GetNScreenOutput();
    double mass      = 0.0;

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

            case ITER: _screenOutputFields[idx_field] = idx_pseudoTime; break;

            case RMS_FLUX:
                _screenOutputFields[idx_field] = blaze::l2Norm( _fluxNew - _flux );
                _flux                          = _fluxNew;
                break;

            default: ErrorMessages::Error( "Screen Output Group not defined!", CURRENT_FUNCTION ); break;
        }
    }
}

void Solver::PrintScreen( std::shared_ptr<spdlog::logger> log ) {
    log->info( "{:03.8f}   {:01.5e} {:01.5e}", _screenOutputFields[0], _screenOutputFields[1], _screenOutputFields[2] );
}

void Solver::Solve() {
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    auto log = spdlog::get( "event" );

    double mass = 0;

    if( rank == 0 ) log->info( "{:10}   {:10}", "t", "mass" );

    // Loop over energies (pseudo-time of continuous slowing down approach)
    for( unsigned idx_energy = 0; idx_energy < _nEnergies; idx_energy++ ) {

        // --- Prepare Boundaries and temp variables
        IterPreprocessing();

        // --- Compute Fluxes ---
        FluxUpdate( _solNew );

        // --- Finite Volume Update ---
        FVMUpdate( _solNew, idx_energy );

        // --- Postprocessing ---
        IterPostprocessing();

        // --- VTK and CSV Output ---
        mass = WriteOutputFields( idx_energy );
        Save( idx_energy );

        // --- Screen Output ---
        WriteScreenOutputFields( idx_energy );
        PrintScreen( log );
    }
}
