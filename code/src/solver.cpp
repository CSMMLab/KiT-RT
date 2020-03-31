#include "solver.h"
#include "snsolver.h"

Solver::Solver( Settings* settings ) : _NCells( 1 ), _NTimeSteps( 1 ) {
    // @TODO save parameters from settings class

    // build quadrature object and store quadrature points and weights
    Quadrature* q = Quadrature::CreateQuadrature( settings->GetQuadName(), settings->GetQuadOrder() );
    _quadPoints   = q->GetPoints();
    _weights      = q->GetWeights();
    _nq           = q->GetNq();

    // setup angular flux array (maybe directly call SetupIC() from physics class? )
    _psi = Matrix( _NCells, _nq, 0.0 );

    // @TODO create object mesh and store _cells, ...

    // @TODO create object quadrature and store _quadPoints, ...

    // @TODO create object quadrature and store _quadPoints, ...

    // load density and stopping power data
    LoadPatientDensity( "test" );
    LoadStoppingPower( "test" );
}

void Solver::LoadPatientDensity( std::string fileName ) {
    // @TODO
}

void Solver::LoadStoppingPower( std::string fileName ) {
    // @TODO
}

Solver* Solver::Create( Settings* settings ) { return new SNSolver( settings ); }
