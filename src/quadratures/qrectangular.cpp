#include "quadratures/qrectangular.hpp"
#include "common/config.hpp"
#include "toolboxes/errormessages.hpp"

#include <iostream>

QRectangular::QRectangular( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    SetNq();
    SetPointsAndWeights();
    _supportedDimensions = { 3 };
}

void QRectangular::SetNq() { _nq = pow( GetOrder(), 3 ); }

void QRectangular::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

void QRectangular::SetPointsAndWeights() {
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

    // Create the grid in karthesian and spherical coordinate
    double global_idx;
    for( unsigned i_x = 0; i_x < _order; i_x++ ) {
        for( unsigned i_y = 0; i_y < _order; i_y++ ) {
            for( unsigned i_z = 0; i_z < _order; i_z++ ) {
                global_idx                  = i_x * _order * _order + i_y * _order + i_z;
                _pointsKarth[global_idx][0] = nodes1D[i_x];
                _pointsKarth[global_idx][1] = nodes1D[i_y];
                _pointsKarth[global_idx][2] = nodes1D[i_z];

                _pointsSphere[global_idx][0] = nodes1D[i_z];    // mu = z due to consistency with 3D case
                if( nodes1D[i_y] >= 0 )
                    _pointsSphere[global_idx][1] = acos( nodes1D[i_x] );    // phi
                else
                    _pointsSphere[global_idx][1] = 2.0 * M_PI - acos( nodes1D[i_x] );    // phi
                _pointsSphere[global_idx][2] =
                    sqrt( nodes1D[i_x] * nodes1D[i_x] + nodes1D[i_y] * nodes1D[i_y] + nodes1D[i_z] * nodes1D[i_z] );    // r radius

                _weights[global_idx] = weights1D[i_x] * weights1D[i_y] * weights1D[i_z];    // Equal weights
            }
        }
    }
}

void QRectangular::ScalePointsAndWeights( double velocityScaling ) {
    // Scale from [-1,1] to [-velocityScaling,velocityScaling] in 1D
    // Scale radius of velocity sphere with velocityScaling in 2D and 3D
    if( !_settings ) {
        ErrorMessages::Error( "This function is only available with an active settings file.", CURRENT_FUNCTION );
    }
    double x, y, z;
    _weights = _weights * velocityScaling;
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        _pointsKarth[idx_quad][0]  = _pointsKarth[idx_quad][0] * velocityScaling;
        _pointsKarth[idx_quad][1]  = _pointsKarth[idx_quad][1] * velocityScaling;
        _pointsKarth[idx_quad][2]  = _pointsKarth[idx_quad][2] * velocityScaling;
        x                          = _pointsKarth[idx_quad][0];
        y                          = _pointsKarth[idx_quad][1];
        z                          = _pointsKarth[idx_quad][2];
        _pointsSphere[idx_quad][0] = z;    // mu = z due to consistency with 3D case
        if( y >= 0 )
            _pointsSphere[idx_quad][1] = acos( x );    // phi
        else
            _pointsSphere[idx_quad][1] = 2.0 * M_PI - acos( x );       // phi
        _pointsSphere[idx_quad][2] = sqrt( x * x + y * y + z * z );    // r radius
    }
}

QRectangular2D::QRectangular2D( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    SetNq();
    CheckOrder();
    SetPointsAndWeights();
    _supportedDimensions = { 2 };
}

void QRectangular2D::SetNq() { _nq = pow( GetOrder(), 2 ); }

void QRectangular2D::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

void QRectangular2D::CheckOrder() {
    if( _order % 2 == 1 ) {    // order needs to be even
        ErrorMessages::Error( "ERROR! Order " + std::to_string( _order ) + " for " + GetName() + " not available. \n Order must be an even number. ",
                              CURRENT_FUNCTION );
    }
}

void QRectangular2D::SetPointsAndWeights() {
    Vector nodes1D( _order, 0.0 ), weights1D( _order, 0.0 );

    // Construct points on x axis
    for( unsigned i = 0; i < _order; ++i ) {
        nodes1D[i]   = -1 + ( i + 0.5 ) * 2 / _order;
        weights1D[i] = 2.0 / double( _order );    // integrate from -1 to 1
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

    // Create the grid in karthesian and spherical coordinates
    for( unsigned i_x = 0; i_x < _order; i_x++ ) {
        for( unsigned i_y = 0; i_y < _order; i_y++ ) {
            _pointsKarth[i_x * _order + i_y][0] = nodes1D[i_x];
            _pointsKarth[i_x * _order + i_y][1] = nodes1D[i_y];
            _pointsKarth[i_x * _order + i_y][2] = 0.0;

            _pointsSphere[i_x * _order + i_y][0] = 0.0;    // mu = z due to consistency with 3D case
            if( nodes1D[i_y] >= 0 )
                _pointsSphere[i_x * _order + i_y][1] = acos( nodes1D[i_x] );    // phi
            else
                _pointsSphere[i_x * _order + i_y][1] = 2.0 * M_PI - acos( nodes1D[i_x] );                                // phi
            _pointsSphere[i_x * _order + i_y][2] = sqrt( nodes1D[i_x] * nodes1D[i_x] + nodes1D[i_y] * nodes1D[i_y] );    // r radius

            _weights[i_x * _order + i_y] = weights1D[i_x] * weights1D[i_y];    // Equal weights
        }
    }
}

void QRectangular2D::ScalePointsAndWeights( double velocityScaling ) {
    // Scale from [-1,1] to [-velocityScaling,velocityScaling] in 1D
    // Scale radius of velocity sphere with velocityScaling in 2D and 3D
    if( !_settings ) {
        ErrorMessages::Error( "This function is only available with an active settings file.", CURRENT_FUNCTION );
    }
    _weights = _weights * velocityScaling * velocityScaling;
    double x, y, z;
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        _pointsKarth[idx_quad][0] = _pointsKarth[idx_quad][0] * velocityScaling;
        _pointsKarth[idx_quad][1] = _pointsKarth[idx_quad][1] * velocityScaling;
        _pointsKarth[idx_quad][2] = _pointsKarth[idx_quad][2] * velocityScaling;
        x                         = _pointsKarth[idx_quad][0];
        y                         = _pointsKarth[idx_quad][1];
        z                         = _pointsKarth[idx_quad][2];

        _pointsSphere[idx_quad][0] = z;                                // my = z
        _pointsSphere[idx_quad][1] = atan2( y, x );                    // phi in [-pi,pi]
        _pointsSphere[idx_quad][2] = sqrt( x * x + y * y + z * z );    // radius r

        // adapt intervall s.t. phi in [0,2pi]
        if( _pointsSphere[idx_quad][1] < 0 ) {
            _pointsSphere[idx_quad][1] = 2 * M_PI + _pointsSphere[idx_quad][1];
        }
    }
}

QRectangular1D::QRectangular1D( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    SetNq();
    SetPointsAndWeights();
    _supportedDimensions = { 1 };
}

void QRectangular1D::SetNq() { _nq = GetOrder(); }

void QRectangular1D::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

void QRectangular1D::SetPointsAndWeights() {
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

void QRectangular1D::ScalePointsAndWeights( double velocityScaling ) {
    // Scale from [-1,1] to [-velocityScaling,velocityScaling] in 1D
    if( !_settings ) {
        ErrorMessages::Error( "This function is only available with an active settings file.", CURRENT_FUNCTION );
    }

    _weights = _weights * velocityScaling;
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        _pointsKarth[idx_quad][0]  = _pointsKarth[idx_quad][0] * velocityScaling;
        _pointsSphere[idx_quad][0] = _pointsSphere[idx_quad][0] * velocityScaling;    // scale radius
    }
}
