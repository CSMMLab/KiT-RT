#include "quadratures/qmidpointtensorized.hpp"
#include "common/config.hpp"
#include "toolboxes/errormessages.hpp"

#include <iostream>

QMidpointTensorized::QMidpointTensorized( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    SetNq();
    SetPointsAndWeights();
    _supportedDimensions = { 3 };
}

void QMidpointTensorized::SetNq() { _nq = 2 * pow( GetOrder(), 2 ); }

void QMidpointTensorized::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

void QMidpointTensorized::SetPointsAndWeights() {
    Vector nodes1D( _order, 0.0 ), weights1D( _order, 0.0 );

    // Construct points on x axis
    for( unsigned i = 0; i < _order; ++i ) {
        nodes1D[i]   = -1 + ( i + 0.5 ) * 2 / _order;
        weights1D[i] = 2.0 / double( _order );    // integrate from -1 to 1
    }

    // setup equidistant angle phi around z axis
    Vector phi( 2 * _order );
    for( unsigned i = 0; i < 2 * _order; ++i ) {
        phi[i] = ( i + 0.5 ) * M_PI / _order;
    }

    // resize points and weights
    _pointsKarth.resize( _nq );
    _pointsSphere.resize( _nq );
    for( auto& p : _pointsKarth ) {
        p.resize( 3 );
    }
    for( auto& p : _pointsSphere ) {
        p.resize( 3 );
    }

    _weights.resize( _nq );

    // transform tensorized (x,y,z)-grid to spherical grid points
    for( unsigned j = 0; j < _order; ++j ) {
        for( unsigned i = 0; i < 2 * _order; ++i ) {
            _pointsKarth[j * ( 2 * _order ) + i][0]  = sqrt( 1 - nodes1D[j] * nodes1D[j] ) * cos( phi[i] );
            _pointsKarth[j * ( 2 * _order ) + i][1]  = sqrt( 1 - nodes1D[j] * nodes1D[j] ) * sin( phi[i] );
            _pointsKarth[j * ( 2 * _order ) + i][2]  = nodes1D[j];
            _pointsSphere[j * ( 2 * _order ) + i][0] = nodes1D[j];    // my
            _pointsSphere[j * ( 2 * _order ) + i][1] = phi[i];        // phi
            _pointsSphere[j * ( 2 * _order ) + i][2] = 1.0;           // radius r
            _weights[j * ( 2 * _order ) + i]         = M_PI / _order * weights1D[j];
        }
    }
}

QMidpointTensorized2D::QMidpointTensorized2D( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    SetNq();
    CheckOrder();
    SetPointsAndWeights();
    _supportedDimensions = { 2 };
}

void QMidpointTensorized2D::SetNq() { _nq = pow( GetOrder(), 2 ); }

void QMidpointTensorized2D::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

void QMidpointTensorized2D::CheckOrder() {
    if( _order % 2 == 1 ) {    // order needs to be even
        ErrorMessages::Error( "ERROR! Order " + std::to_string( _order ) + " for " + GetName() + " not available. \n Order must be an even number. ",
                              CURRENT_FUNCTION );
    }
}

void QMidpointTensorized2D::SetPointsAndWeights() {
    Vector nodes1D( _order, 0.0 ), weights1D( _order, 0.0 );

    // Construct points on x axis
    for( unsigned i = 0; i < _order; ++i ) {
        nodes1D[i]   = -1 + ( i + 0.5 ) * 2 / _order;
        weights1D[i] = 2.0 / double( _order );    // integrate from -1 to 1
    }

    // setup equidistant angle phi around z axis
    Vector phi( 2 * _order );
    for( unsigned i = 0; i < 2 * _order; ++i ) {
        phi[i] = ( i + 0.5 ) * M_PI / _order;
    }

    // resize points and weights
    _pointsKarth.resize( _nq );
    _pointsSphere.resize( _nq );
    for( auto& p : _pointsKarth ) {
        p.resize( 3 );
    }
    for( auto& p : _pointsSphere ) {
        p.resize( 3 );
    }

    _weights.resize( _nq );
    unsigned range = std::floor( _order / 2.0 );    // comment (steffen): Only half of the points, due to projection

    // transform tensorized (x,y,z)-grid to spherical grid points
    for( unsigned j = 0; j < range; ++j ) {
        for( unsigned i = 0; i < 2 * _order; ++i ) {
            _pointsKarth[j * ( 2 * _order ) + i][0]  = nodes1D[j] * cos( phi[i] );
            _pointsKarth[j * ( 2 * _order ) + i][1]  = nodes1D[j] * sin( phi[i] );
            _pointsKarth[j * ( 2 * _order ) + i][2]  = 0.0;
            _pointsSphere[j * ( 2 * _order ) + i][0] = _pointsKarth[j * ( 2 * _order ) + i][2];    // my = z
            _pointsSphere[j * ( 2 * _order ) + i][1] =
                atan2( _pointsKarth[j * ( 2 * _order ) + i][1], _pointsKarth[j * ( 2 * _order ) + i][0] );    // phi in [-pi,pi]
            _pointsSphere[j * ( 2 * _order ) + i][2] =
                sqrt( _pointsKarth[j * ( 2 * _order ) + i][0] * _pointsKarth[j * ( 2 * _order ) + i][0] +
                      _pointsKarth[j * ( 2 * _order ) + i][1] * _pointsKarth[j * ( 2 * _order ) + i][1] );    // radius r

            // adapt intervall s.t. phi in [0,2pi]
            if( _pointsSphere[j * ( 2 * _order ) + i][1] < 0 ) {
                _pointsSphere[j * ( 2 * _order ) + i][1] = 2 * M_PI + _pointsSphere[j * ( 2 * _order ) + i][1];
            }
            _weights[j * ( 2 * _order ) + i] = M_PI / ( 2.0 * _order ) * weights1D[j];
        }
    }
}

QMidpoint1D::QMidpoint1D( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    SetNq();
    SetPointsAndWeights();
    _supportedDimensions = { 1 };
}

void QMidpoint1D::SetNq() { _nq = GetOrder(); }

void QMidpoint1D::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

void QMidpoint1D::SetPointsAndWeights() {
    Vector nodes1D( _order, 0.0 ), weights1D( _order, 0.0 );

    // Construct points on x axis
    for( unsigned i = 0; i < _order; ++i ) {
        nodes1D[i]   = -1 + ( i + 0.5 ) * 2.0 / (double)_order;
        weights1D[i] = 2.0 / double( _order );    // integrate from -1 to 1
    }

    // setup equidistant angle phi around z axis
    Vector phi( 2 * _order );
    for( unsigned i = 0; i < 2 * _order; ++i ) {
        phi[i] = ( i + 0.5 ) * M_PI / double( _order );
    }

    // resize points and weights
    _pointsKarth.resize( _nq );
    _weights.resize( _nq );
    unsigned dim = 1;
    for( unsigned k = 0; k < _nq; ++k ) {
        _pointsKarth[k].resize( dim );
        _pointsKarth[k][0] = nodes1D[k];
        _weights[k]        = weights1D[k];
    }
    _pointsSphere = _pointsKarth;
}
