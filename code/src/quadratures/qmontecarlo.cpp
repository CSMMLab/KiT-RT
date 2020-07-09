#include "../../include/quadratures/qmontecarlo.h"

QMonteCarlo::QMonteCarlo( unsigned order ) : QuadratureBase( order ) {
    SetName();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
}

void QMonteCarlo::SetPointsAndWeights() {
    // Nq random points on the sphere.
    _points.resize( GetNq() );
    _pointsSphere.resize( GetNq() );
    std::default_random_engine generator;
    std::normal_distribution<double> distribution( 0.0, 1.0 );
    std::uniform_real_distribution<double> distributionUni( 0.0, 1.0 );
    // Tranlation in spherical coordinates
    double my, phi, x, y, z, temp, norm;

    for( unsigned i = 0; i < GetNq(); i++ ) {
        // https://mathworld.wolfram.com/SpherePointPicking.html describes how to generate random points on the sphere.
        x    = distribution( generator );
        y    = distribution( generator );
        z    = distribution( generator );
        norm = sqrt( x * x + y * y + z * z );
        Vector point( { x / norm, y / norm, z / norm } );
        _points[i] = point;

        // Version for spherical coordinates (not 1:1 correspondence to karthesian in terms of distribution!)
        temp = distributionUni( generator );
        my   = 2 * temp - 1;
        phi  = 2 * M_PI * temp;

        Vector pointSphere( { my, phi } );
        _pointsSphere[i] = pointSphere;
    }

    // Equal weights
    _weights = Vector( GetNq(), 4.0 * M_PI / GetNq() );
}

VectorVector QMonteCarlo::GetPointsSphere() const { return _pointsSphere; }

void QMonteCarlo::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}
