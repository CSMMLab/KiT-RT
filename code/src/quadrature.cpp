#include "quadrature.h"
#include "qmontecarlo.h"
#include "qgausslegendretensorized.h"
#include "qlevelsymmetric.h"
#include "qldfesa.h"
#include "qlebedev.h"

Quadrature::Quadrature( unsigned order ) : _order( order ) {}

Quadrature* Quadrature::CreateQuadrature( QUAD_NAME name, unsigned order ) {

    switch (name){
        case QUAD_MonteCarlo:               return new QMonteCarlo( order );
        case QUAD_GaussLegendreTensorized:  return new QGaussLegendreTensorized( order );
        case QUAD_LevelSymmetric:           return new QLevelSymmetric( order );
        case QUAD_LDFESA:                   return new QLDFESA(order);
        case QUAD_Lebedev:                  return new QLebedev(order);
        default:                            return new QMonteCarlo( order ); // Use MonteCarlo as dummy
    }
}

double Quadrature::Integrate( double( f )( double x0, double x1, double x2 ) ) {
    double result = 0;
    for( unsigned i = 0; i < _nq; i++ ) {
        double x = _points[i][0];
        double y = _points[i][1];
        double z = _points[i][2];
        double w = _weights[i];
        result += w * f( x, y, z );
    }
    return result;
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
