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
    //_s       = _problem->GetStoppingPower( _energies );
    _sigmaT = _problem->GetTotalXS( _energies );
    _sigmaS = _problem->GetScatteringXS( _energies );
    _Q      = _problem->GetExternalSource( _energies );

    // setup numerical flux
    _g = NumericalFlux::Create();

    // boundary type
    _boundaryCells = _mesh->GetBoundaryTypes();

    // Solver Output
    _solverOutput.resize( _nCells );    // Currently only Flux
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

void Solver::Save() const {
    std::vector<std::string> fieldNames{ "flux" };
    std::vector<double> flux;
    flux.resize( _nCells );

    for( unsigned i = 0; i < _nCells; ++i ) {
        flux[i] = _sol[i][0];
    }
    std::vector<std::vector<double>> scalarField( 1, flux );
    std::vector<std::vector<std::vector<double>>> results{ scalarField };
    ExportVTK( _settings->GetOutputFile() + "_" + std::to_string( _nEnergies ), results, fieldNames, _mesh );
}
