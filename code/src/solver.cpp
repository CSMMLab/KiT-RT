#include "solver.h"
#include "../../include/quadratures/quadraturebase.h"
#include "io.h"
#include "mesh.h"
#include "snsolver.h"

Solver::Solver( Config* settings ) : _settings( settings ) {
    // @TODO save parameters from settings class

    // build mesh and store all relevant information
    _mesh      = LoadSU2MeshFromFile( settings );
    _areas     = _mesh->GetCellAreas();
    _neighbors = _mesh->GetNeighbours();
    _normals   = _mesh->GetNormals();
    _nCells    = _mesh->GetNumCells();
    _settings->SetNCells( _nCells );

    // build quadrature object and store quadrature points and weights
    QuadratureBase* q = QuadratureBase::CreateQuadrature( settings->GetQuadName(), settings->GetQuadOrder() );
    _quadPoints       = q->GetPoints();
    _weights          = q->GetWeights();
    _nq               = q->GetNq();
    _settings->SetNQuadPoints( _nq );
    ScatteringKernel* k = ScatteringKernel::CreateScatteringKernel( settings->GetKernelName(), _nq );
    _scatteringKernel   = k->GetScatteringKernel();

    // set time step
    _dE        = ComputeTimeStep( settings->GetCFL() );
    _nEnergies = unsigned( settings->GetTEnd() / _dE );
    for( unsigned i = 0; i < _nEnergies; ++i ) _energies.push_back( ( i + 1 ) * _dE );

    // setup problem
    _problem = ProblemBase::Create( _settings, _mesh );
    _psi     = _problem->SetupIC();
    _s       = _problem->GetStoppingPower( _energies );
    _sigmaT  = _problem->GetTotalXS( _energies );
    _sigmaS  = _problem->GetScatteringXS( _energies );
    _Q       = _problem->GetExternalSource( _energies );

    // setup numerical flux
    _g = NumericalFlux::Create( settings );

    // boundary type
    _boundaryCells = _mesh->GetBoundaryTypes();
}

double Solver::ComputeTimeStep( double cfl ) const {
    double maxEdge = -1.0;
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned l = 0; l < _normals[j].size(); ++l ) {
            double currentEdge = _areas[j] / norm( _normals[j][l] );
            if( currentEdge > maxEdge ) maxEdge = currentEdge;
        }
    }
    return cfl * maxEdge;
}

Solver* Solver::Create( Config* settings ) { return new SNSolver( settings ); }
