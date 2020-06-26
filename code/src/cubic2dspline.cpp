#include "cubic2dspline.h"

#include <iostream>

Cubic2DSpline::Cubic2DSpline( const Vector& x, const Vector& y, const Matrix& data ) : _x( x ), _y( y ), _data( data ) {}
Cubic2DSpline::~Cubic2DSpline() {}

inline double Cubic2DSpline::interpolate1D( double param[4], double x ) {
    return ( param[1] + 0.5 * x *
                            ( param[2] - param[0] +
                              x * ( 2.0 * param[0] - 5.0 * param[1] + 4.0 * param[2] - param[3] +
                                    x * ( 3.0 * ( param[1] - param[2] ) + param[3] - param[0] ) ) ) );
}

double Cubic2DSpline::operator()( double x, double y ) {
    unsigned xId = indexOfClosestValue( x, _x );
    unsigned yId = indexOfClosestValue( y, _y );

    // store all 16 interpolation points needed
    double points[4][4];
    for( int i = -1; i < 3; ++i ) {
        unsigned idx_y;
        idx_y = yId + i < 0 ? 0 : yId + i;
        idx_y = yId + i > _data.rows() - 1 ? _data.rows() - 1 : yId + i;
        for( int j = -1; j < 3; ++j ) {
            unsigned idx_x;
            idx_x = xId + j < 0 ? 0 : xId + j;
            idx_x = xId + j > _data.columns() - 1 ? _data.columns() - 1 : xId + j;

            points[i + 1][j + 1] = _data( idx_x, idx_y );
        }
    }

    // rescale data to [0,1]
    double t = ( x - _x[xId] ) / ( _x[xId + 1] - _x[xId] );
    double u = ( y - _y[yId] ) / ( _y[yId + 1] - _y[yId] );

    // first interpolate in x-direction and store the results on which the final interpolation in y will be done
    double interpolationBuffer[4];
    interpolationBuffer[0] = interpolate1D( points[0], t );
    interpolationBuffer[1] = interpolate1D( points[1], t );
    interpolationBuffer[2] = interpolate1D( points[2], t );
    interpolationBuffer[3] = interpolate1D( points[3], t );
    return interpolate1D( interpolationBuffer, u );
}

unsigned Cubic2DSpline::indexOfClosestValue( double value, const Vector& v ) {
    auto i = std::min_element( begin( v ), end( v ), [=]( double x, double y ) { return abs( x - value ) < abs( y - value ); } );
    return std::distance( begin( v ), i );
}
