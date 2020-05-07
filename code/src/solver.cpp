#include "solver.h"
#include "../../include/quadratures/quadrature.h"
#include "io.h"
#include "mesh.h"
#include "snsolver.h"

Solver::Solver( Settings* settings ) : _settings( settings ) {
    // @TODO save parameters from settings class

    // std::cout << "In Solver..." << std::endl;

    // build quadrature object and store quadrature points and weights
    Quadrature* q = Quadrature::CreateQuadrature( settings->GetQuadName(), settings->GetQuadOrder() );
    _quadPoints   = q->GetPoints();
    _weights      = q->GetWeights();
    _nq           = q->GetNq();

    // build mesh and store all relevant information

    _mesh      = LoadSU2MeshFromFile( settings );
    _areas     = _mesh->GetCellAreas();
    _neighbors = _mesh->GetNeighbours();
    _normals   = _mesh->GetNormals();
    _nCells    = _mesh->GetNumCells();

    // std::cout << "After Mesh..." << std::endl;

    // setup angular flux array (maybe directly call SetupIC() from physics class? )
    _psi = std::vector( _nCells, Vector( _nq, 1e-7 ) );

    Vector midPoint( 2, 0.6 );
    auto nodes = _mesh->GetCellMidPoints();
    // write some IC
    for( unsigned j = 0; j < nodes.size(); ++j ) {

        // std::cout << norm( nodes[j] - midPoint ) << std::endl;
        if( norm( nodes[j] - midPoint ) <= 0.1 ) {
            // std::cout << nodes[j] << std::endl;
            for( unsigned k = 0; k < _nq; ++k ) {
                _psi[j][k] = 1e-7;
            }
        }
    }

    // set time step
    _dt         = ComputeTimeStep( settings->GetCFL() );
    _nTimeSteps = unsigned( settings->GetTEnd() / _dt );

    // TODO: load density and stopping power data
    LoadPatientDensity( "test" );
    LoadStoppingPower( "test" );
    LoadSigmaS( "test" );
    LoadSigmaT( "test" );

    // TODO: setup initial condtion (distribution at initial energy for CSD)
    SetupIC();

    // setup numerical flux
    _g = NumericalFlux::Create( settings );

    // boundary type
    _boundaryCells = _mesh->GetBoundaryCellArray();
}

double Solver::ComputeTimeStep( double cfl ) const {
    double maxEdge = -1.0;
    // std::cout << _nCells << std::endl;
    // std::cout << _areas.size() << std::endl;
    // std::cout << _normals.size() << std::endl;
    for( unsigned j = 0; j < _nCells; ++j ) {
        // std::cout << j;
        // std::cout << " " << _areas[j] << " " << _normals[j].size() << std::endl;
        for( unsigned l = 0; l < _normals[j].size(); ++l ) {
            // std::cout << _normals[j][l] << std::endl;
            double currentEdge = _areas[j] / norm( _normals[j][l] );
            if( currentEdge > maxEdge ) maxEdge = currentEdge;
        }
    }
    return cfl * maxEdge;
}

void Solver::LoadPatientDensity( std::string fileName ) {
    // @TODO
    // _density = ... -> dim(_density) = _nCells
}

void Solver::LoadStoppingPower( std::string fileName ) {
    // @TODO
    //_sH20 = ... -> dim(_sH20) = _nTimeSteps
}

void Solver::LoadSigmaS( std::string fileName ) {
    // @TODO
    //_sigmaSH20 = ... -> dim(_sigmaSH20) = (_nTimeSteps,_nq,_nq)
    _sigmaSH20 = std::vector( _nTimeSteps, Matrix( _nq, _nq, 0.0 ) );    // TODO: double check this!
}

void Solver::LoadSigmaT( std::string fileName ) {
    // @TODO
    //_sigmaTH20 = ... -> dim(_sigmaTH20) = _nTimeSteps
    _sigmaTH20 = std::vector( _nTimeSteps, 0.0 );
}

void Solver::SetupIC() {
    // @TODO
    // write initial condition
    // _psi = ...
}

Solver* Solver::Create( Settings* settings ) { return new SNSolver( settings ); }
