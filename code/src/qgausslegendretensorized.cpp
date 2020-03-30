#include "qgausslegendretensorized.h"

QGaussLegendreTensorized::QGaussLegendreTensorized( int order ) : Quadrature( order ) {
    SetName( ComputeName() );
    SetNq( ComputeNq() );
    SetPoints( ComputePoints() );
    SetWeights( ComputeWeights() );
    SetConnectivity( ComputeConnectivity() );
}

std::string QGaussLegendreTensorized::ComputeName() { return "Tensorized Gauss-Legendre quadrature."; };

int QGaussLegendreTensorized::ComputeNq() { return pow( GetOrder(), 2 ); };

blaze::DynamicVector<blaze::DynamicVector<double>> QGaussLegendreTensorized::ComputePoints() {
    // Nq random points on the sphere.
    int nq = GetNq();
    blaze::DynamicVector<blaze::DynamicVector<double>> points( nq );

    std::default_random_engine generator;
    std::normal_distribution<double> distribution( 0.0, 1.0 );

    for( int i = 0; i < nq; i++ ) {
        // https://mathworld.wolfram.com/SpherePointPicking.html describes how to generate random points on the sphere.
        double x    = distribution( generator );
        double y    = distribution( generator );
        double z    = distribution( generator );
        double norm = sqrt( x * x + y * y + z * z );
        blaze::DynamicVector<double> point( {x / norm, y / norm, z / norm} );
        points[i] = point;
    }
    return points;
};

blaze::DynamicVector<double> QGaussLegendreTensorized::ComputeWeights() {
    // Equal weights
    int nq = GetNq();
    blaze::DynamicVector<double> weights( nq, 4.0 * M_PI / nq );
    return weights;
};

blaze::DynamicVector<blaze::DynamicVector<int>> QGaussLegendreTensorized::ComputeConnectivity() {
    // Not initialized for this quadrature.
    blaze::DynamicVector<blaze::DynamicVector<int>> connectivity;
    return connectivity;
};
