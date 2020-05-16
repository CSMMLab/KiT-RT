#include "../../include/quadratures/qmontecarlo.h"

QMonteCarlo::QMonteCarlo( unsigned order ) : QuadratureBase( order ) {
    SetName();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
}

void QMonteCarlo::SetPointsAndWeights() {
    // Nq random points on the sphere.
    _points.resize(GetNq());

    std::default_random_engine generator;
    std::normal_distribution<double> distribution( 0.0, 1.0 );

    for( unsigned i = 0; i < GetNq(); i++ ) {
        // https://mathworld.wolfram.com/SpherePointPicking.html describes how to generate random points on the sphere.
        double x    = distribution( generator );
        double y    = distribution( generator );
        double z    = distribution( generator );
        double norm = sqrt( x * x + y * y + z * z );
        Vector point( {x / norm, y / norm, z / norm} );
        _points[i] = point;
    }

    // Equal weights
    _weights = Vector( GetNq(), 4.0 * M_PI / GetNq() );

}

void QMonteCarlo::SetConnectivity() { //TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}
