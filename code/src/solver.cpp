#include "solver.h"
#include "snsolver.h"

Solver::Solver( Settings* settings ) : _Q( 1 ), _NCells( 1 ), _NTimeSteps( 1 ) {
    // @TODO save parameters from settings class

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
