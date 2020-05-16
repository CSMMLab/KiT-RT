#include "quadratures/quadraturebase.h"
#include "quadratures/qmontecarlo.h"
#include "quadratures/qgausslegendretensorized.h"
#include "quadratures/qlevelsymmetric.h"
#include "quadratures/qldfesa.h"
#include "quadratures/qlebedev.h"

QuadratureBase::QuadratureBase( unsigned order ) : _order( order ) {}

QuadratureBase* QuadratureBase::CreateQuadrature( QUAD_NAME name, unsigned order ) {

    switch (name){
        case QUAD_MonteCarlo:               return new QMonteCarlo( order );
        case QUAD_GaussLegendreTensorized:  return new QGaussLegendreTensorized( order );
        case QUAD_LevelSymmetric:           return new QLevelSymmetric( order );
        case QUAD_LDFESA:                   return new QLDFESA(order);
        case QUAD_Lebedev:                  return new QLebedev(order);
        default:                            return new QMonteCarlo( order ); // Use MonteCarlo as dummy
    }
}

double QuadratureBase::Integrate( double( f )( double x0, double x1, double x2 ) ) {
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

double QuadratureBase::SumUpWeights() { return sum( _weights ); }

void QuadratureBase::PrintWeights() {
    for( unsigned i = 0; i < _nq; i++ ) {
        double w = _weights[i];
        std::cout << w << std::endl;
    }
}

void QuadratureBase::PrintPoints() {
    for( unsigned i = 0; i < _nq; i++ ) {
        double x = _points[i][0];
        double y = _points[i][1];
        double z = _points[i][2];
        std::cout << x << ", " << y << ", " << z << std::endl;
    }
}
void QuadratureBase::PrintPointsAndWeights() {
    for( unsigned i = 0; i < _nq; i++ ) {
        double x = _points[i][0];
        double y = _points[i][1];
        double z = _points[i][2];
        double w = _weights[i];
        std::cout << x << ", " << y << ", " << z << ", " << w << std::endl;
    }
}
