#include "quadrature.h"
#include "qmontecarlo.h"

Quadrature::Quadrature( int order ) : _order( order ) {}

Quadrature* Quadrature::CreateQuadrature( std::string name, int order ) {

    if( name == "montecarlo" ) {
        return new QMonteCarlo( order );
    }

    // If nothing has been picked, take this as dummy:
    return new QMonteCarlo( order );
}

double Quadrature::SumUpWeights() { return sum( _weights ); }

void Quadrature::PrintWeights() {
    for( int i = 0; i < _nq; i++ ) {
        std::cout << _weights[i] << std::endl;
    }
}

void Quadrature::PrintPoints() {
    for( int i = 0; i < _nq; i++ ) {
        std::cout << _points[i][0] << ", " << _points[i][1] << ", " << _points[i][2] << std::endl;
    }
}
