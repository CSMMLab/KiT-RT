#include "cubic2dspline.h"

#include <iostream>

Cubic2DSpline::Cubic2DSpline( const Vector& x, const Vector& y, const Matrix& data ) : _x( x ), _y( y ), _data( data ) {
    //_data = addGhostLayers( data );
}
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
    /*
    if( interpolate1D( interpolationBuffer, y ) < 0 || interpolate1D( interpolationBuffer, y ) > 1 ) {
        std::cerr << points[0][0] << " " << points[0][1] << " " << points[0][2] << " " << points[0][3] << std::endl;
        std::cerr << points[1][0] << " " << points[1][1] << " " << points[1][2] << " " << points[1][3] << std::endl;
        std::cerr << points[2][0] << " " << points[2][1] << " " << points[2][2] << " " << points[2][3] << std::endl;
        std::cerr << points[3][0] << " " << points[3][1] << " " << points[3][2] << " " << points[3][3] << std::endl;
        std::cerr << std::endl;
        std::cerr << interpolationBuffer[0] << " " << interpolationBuffer[1] << " " << interpolationBuffer[2] << " " << interpolationBuffer[3] << " "
                  << x << std::endl;
        std::cerr << std::endl << std::endl;
    }
    */
    return std::clamp( interpolate1D( interpolationBuffer, u ), 0.0, 1.0 );
}

Matrix Cubic2DSpline::addGhostLayers( const Matrix& m ) {
    Matrix paddedMatrix( m.rows() + 4, m.columns() + 4, 0.0 );
    auto subMat = submatrix( paddedMatrix, 2u, 2u, m.rows(), m.columns() );
    subMat      = m;

    // x-direction
    for( unsigned i = 0; i < m.rows(); ++i ) {
        paddedMatrix( i, 0 )               = m( i, 0 );
        paddedMatrix( i, 1 )               = m( i, 0 );
        paddedMatrix( i, m.columns() )     = m( i, m.columns() - 1 );
        paddedMatrix( i, m.columns() + 1 ) = m( i, m.columns() - 1 );
    }

    // y-direction
    for( unsigned i = 0; i < m.columns(); ++i ) {
        paddedMatrix( 0, i )            = m( 0, i );
        paddedMatrix( 1, i )            = m( 0, i );
        paddedMatrix( m.rows(), i )     = m( m.rows() - 1, i );
        paddedMatrix( m.rows() + 1, i ) = m( m.rows() - 1, i );
    }

    // ascending diagonal
    paddedMatrix( 0, 0 )                          = m( 0, 0 );
    paddedMatrix( 1, 1 )                          = m( 0, 0 );
    paddedMatrix( m.rows(), m.columns() )         = m( m.rows() - 1, m.columns() - 1 );
    paddedMatrix( m.rows() + 1, m.columns() + 1 ) = m( m.rows() - 1, m.columns() - 1 );

    // decending diagonal
    paddedMatrix( 0, m.columns() + 1 ) = m( 0, m.columns() - 1 );
    paddedMatrix( 1, m.columns() )     = m( 0, m.columns() - 1 );
    paddedMatrix( m.rows(), 0 )        = m( m.rows() - 1, 0 );
    paddedMatrix( m.rows() + 1, 1 )    = m( m.rows() - 1, 0 );

    // untreated off-diagonal terms
    unsigned i = paddedMatrix.rows() - 1;
    unsigned j = paddedMatrix.columns() - 1;

    paddedMatrix( 0, 1 )     = paddedMatrix( 0, 2 );
    paddedMatrix( 1, 0 )     = paddedMatrix( 2, 0 );
    paddedMatrix( i - 1, j ) = paddedMatrix( i - 2, j );
    paddedMatrix( i, j - 1 ) = paddedMatrix( i, j - 2 );
    paddedMatrix( 0, j - 1 ) = paddedMatrix( 0, j - 2 );
    paddedMatrix( 1, j )     = paddedMatrix( 2, j );
    paddedMatrix( i - 1, 0 ) = paddedMatrix( i - 2, 0 );
    paddedMatrix( i, 1 )     = paddedMatrix( i, 2 );

    return paddedMatrix;
}

unsigned Cubic2DSpline::indexOfClosestValue( double value, const Vector& v ) {
    auto i = std::min_element( begin( v ), end( v ), [=]( double x, double y ) { return abs( x - value ) < abs( y - value ); } );
    return std::distance( begin( v ), i );
}
