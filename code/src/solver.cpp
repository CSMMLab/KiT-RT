#include "solver.h"
#include "mesh.h"
#include "../../include/quadratures/quadrature.h"
#include "snsolver.h"

Solver::Solver( Settings* settings ) : _NCells( 1 ) {
    // @TODO save parameters from settings class

    // build quadrature object and store quadrature points and weights
    Quadrature* q = Quadrature::CreateQuadrature( settings->GetQuadName(), settings->GetQuadOrder() );
    _quadPoints   = q->GetPoints();
    _weights      = q->GetWeights();
    _nq           = q->GetNq();

    // setup angular flux array (maybe directly call SetupIC() from physics class? )
    _psi = std::vector( _NCells, Vector( _nq, 0.0 ) );

    // build mesh and store all relevant information
    // Mesh* mesh = Mesh::Create();
    //_areas = mesh->GetAreas();
    //_neighbors = mesh->GetNeighbors();
    //_normals = mesh->GetNormals();

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
}

double Solver::ComputeTimeStep( double cfl ) const {
    double maxEdge = -1.0;
    for( unsigned j = 0; j < _NCells; ++j ) {
        for( unsigned l = 0; l < _normals[j].size(); ++j ) {
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
}

void Solver::LoadSigmaT( std::string fileName ) {
    // @TODO
    //_sigmaTH20 = ... -> dim(_sigmaTH20) = _nTimeSteps
}

void Solver::SetupIC() {
    // @TODO
    // write initial condition
    // _psi = ...
}

Solver* Solver::Create( Settings* settings ) { return new SNSolver( settings ); }
