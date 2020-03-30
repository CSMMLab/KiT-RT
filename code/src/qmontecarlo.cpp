#include "qmontecarlo.h"

QMonteCarlo::QMonteCarlo( unsigned order ) : Quadrature( order ) {
    SetName( ComputeName() );
    SetNq( ComputeNq() );
    SetPoints( ComputePoints() );
    SetWeights( ComputeWeights() );
    SetConnectivity( ComputeConnectivity() );
};

std::string QMonteCarlo::ComputeName() {
    // Human readable description.
    return "Monte Carlo Quadrature.";
};

unsigned QMonteCarlo::ComputeNq() {
    // For Monte Carlo Quadrature, nq = order.
    return GetOrder();
};

VectorVector QMonteCarlo::ComputePoints() {
    // Nq random points on the sphere.
    unsigned nq = GetNq();
    VectorVector points( nq );

    std::default_random_engine generator;
    std::normal_distribution<double> distribution( 0.0, 1.0 );

    for( unsigned i = 0; i < nq; i++ ) {
        // https://mathworld.wolfram.com/SpherePointPicking.html describes how to generate random points on the sphere.
        double x    = distribution( generator );
        double y    = distribution( generator );
        double z    = distribution( generator );
        double norm = sqrt( x * x + y * y + z * z );
        Vector point( {x / norm, y / norm, z / norm} );
        points[i] = point;
    }
    return points;
};

Vector QMonteCarlo::ComputeWeights() {
    // Equal weights
    unsigned nq = GetNq();
    Vector weights( nq, 4.0 * M_PI / nq );
    return weights;
};

VectorVectorU QMonteCarlo::ComputeConnectivity() {
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    return connectivity;
};
