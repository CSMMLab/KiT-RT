#include "quadratures/qmontecarlo.hpp"

#include <math.h>

QMonteCarlo::QMonteCarlo( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
    _supportedDimensions = { 1, 2, 3 };
}
QMonteCarlo::QMonteCarlo( unsigned quadOrder ) : QuadratureBase( quadOrder ) {
    SetName();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
    _supportedDimensions = { 1, 2, 3 };
}

void QMonteCarlo::SetPointsAndWeights() {
    // Nq random points on the sphere.
    _pointsKarth.resize( GetNq(), Vector( 3, 0.0 ) );
    _pointsSphere.resize( GetNq(), Vector( 2, 0.0 ) );
    std::default_random_engine generator;
    std::normal_distribution<double> distribution( 0.0, 1.0 );
    std::uniform_real_distribution<double> distributionUni( 0.0, 1.0 );
    // Tranlation in spherical coordinates
    double x, y, z, norm;

    for( unsigned i = 0; i < GetNq(); i++ ) {
        // https://mathworld.wolfram.com/SpherePointPicking.html describes how to generate random points on the sphere.
        x             = distribution( generator );
        y             = distribution( generator );
        z             = distribution( generator );
        norm          = sqrt( x * x + y * y + z * z );
        _pointsKarth[i][0] = x / norm;
        _pointsKarth[i][1] = y / norm;
        _pointsKarth[i][2] = z / norm;

        // Version for spherical coordinates (not 1:1 correspondence to karthesian in terms of distribution!)

        // temp = distributionUni( generator );
        // my   = 2 * temp - 1;
        // phi  = 2 * M_PI * temp;
        //
        // Vector pointSphere( { my, phi } );
        // _pointsSphere[i] = pointSphere;
    }
    // Transform _points to _pointsSphere ==>transform (x,y,z) into (my,phi)
    for( unsigned idx = 0; idx < _nq; idx++ ) {
        _pointsSphere[idx].resize( 2 );                                       // (my,phi)
        _pointsSphere[idx][0] = _pointsKarth[idx][2];                              // my = z
        _pointsSphere[idx][1] = atan2( _pointsKarth[idx][1], _pointsKarth[idx][0] );    // phi in [-pi,pi]

        // adapt intervall s.t. phi in [0,2pi]
        if( _pointsSphere[idx][1] < 0 ) {
            _pointsSphere[idx][1] = 2 * M_PI + _pointsSphere[idx][1];
        }
    }

    // Equal weights
    _weights = Vector( _nq, 4.0 * M_PI / ( (double)_nq ) );
}

void QMonteCarlo::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}
