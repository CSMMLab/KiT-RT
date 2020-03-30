#include "qmontecarlo.h"

QMonteCarlo::QMonteCarlo( int order ) : Quadrature( order ) {
    SetName( ComputeName() );
    SetNq( ComputeNq() );
    SetPoints( ComputePoints() );
    SetWeights( ComputeWeights() );
    SetConnectivity( ComputeConnectivity() );
};

std::string QMonteCarlo::ComputeName() {
    // Human readable description.
    return "Monte Carlo Quadrature";
};

int QMonteCarlo::ComputeNq() {
    // For Monte Carlo Quadrature, nq = order.
    return GetOrder();
};

blaze::DynamicVector<blaze::DynamicVector<double>> QMonteCarlo::ComputePoints() {
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

blaze::DynamicVector<double> QMonteCarlo::ComputeWeights() {
    // Equal weights
    int nq = GetNq();
    blaze::DynamicVector<double> weights( nq, 4.0 * M_PI / nq );
    return weights;
};

blaze::DynamicVector<blaze::DynamicVector<int>> QMonteCarlo::ComputeConnectivity() {
    // Not initialized for this quadrature.
    blaze::DynamicVector<blaze::DynamicVector<int>> connectivity;
    return connectivity;
};
