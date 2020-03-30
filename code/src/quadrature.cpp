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

double Quadrature::SumUpWeights() {
    double w = 0;
    for( int i = 0; i < _nq; i++ ) {
        w += _weights[i];
    }
    return w;
}
void Quadrature::PrintWeights() {
    for( int i = 0; i < _nq; i++ ) {
        std::cout << _weights[i] << std::endl;
    }
}
