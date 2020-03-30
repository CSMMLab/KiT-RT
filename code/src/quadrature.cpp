#include "quadrature.h"
#include "qmontecarlo.h"

Quadrature::Quadrature( int order ) : _order( order ) {}

Quadrature* GetQuadrature( std::string name, int order ) {
    if( name == "montecarlo" ) {
        return new QMonteCarlo( order );
    }
    return new QMonteCarlo( order );
}
