#include "qgausslegendretensorized.h"

QGaussLegendreTensorized::QGaussLegendreTensorized( unsigned order ) : Quadrature( order ) {
    SetName( ComputeName() );
    SetNq( ComputeNq() );
    SetPoints( ComputePoints() );
    SetWeights( ComputeWeights() );
    SetConnectivity( ComputeConnectivity() );
}

std::string QGaussLegendreTensorized::ComputeName() { return "Tensorized Gauss-Legendre quadrature."; };

unsigned QGaussLegendreTensorized::ComputeNq() { return pow( GetOrder(), 2 ); };

VectorVector QGaussLegendreTensorized::ComputePoints() {
    // Nq random points on the sphere.
    unsigned nq = GetNq();
    VectorVector points( nq );
    return points;
};

Vector QGaussLegendreTensorized::ComputeWeights() {
    // Equal weights
    unsigned nq = GetNq();
    Vector weights( nq, 4.0 * M_PI / nq );
    return weights;
};

VectorVectorU QGaussLegendreTensorized::ComputeConnectivity() {
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    return connectivity;
};
