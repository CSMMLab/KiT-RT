#include "../../include/quadratures/qgausslegendretensorized.h"

QGaussLegendreTensorized::QGaussLegendreTensorized( unsigned order ) : QuadratureBase( order ) {
    SetName();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
}

void QGaussLegendreTensorized::SetPointsAndWeights() {    // TODO
    // Compute Points
    // Nq random points on the sphere.
    _points = VectorVector( GetNq() );

    // Compute Weights
    _weights = Vector( GetNq(), 4.0 * M_PI / GetNq() );
}

void QGaussLegendreTensorized::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}
