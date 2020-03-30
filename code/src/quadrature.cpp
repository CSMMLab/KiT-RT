#include "quadrature.h"
#include "qmontecarlo.h"

Quadrature::Quadrature( unsigned order ) : _order( order ) {}

Quadrature* Quadrature::CreateQuadrature( std::string name, unsigned order ) {

    if( name == "montecarlo" ) {
        return new QMonteCarlo( order );
    }

    // If nothing has been picked, take this as dummy:
    return new QMonteCarlo( order );
}

double Quadrature::SumUpWeights() { return sum( _weights ); }

void Quadrature::PrintWeights() {
    for( unsigned i = 0; i < _nq; i++ ) {
        double w = _weights[i];
        std::cout << w << std::endl;
    }
}

void Quadrature::PrintPoints() {
    for( unsigned i = 0; i < _nq; i++ ) {
        double x = _points[i][0];
        double y = _points[i][1];
        double z = _points[i][2];
        std::cout << x << ", " << y << ", " << z << std::endl;
    }
}
void Quadrature::PrintPointsAndWeights() {
    for( unsigned i = 0; i < _nq; i++ ) {
        double x = _points[i][0];
        double y = _points[i][1];
        double z = _points[i][2];
        double w = _weights[i];
        std::cout << x << ", " << y << ", " << z << ", " << w << std::endl;
    }
}
