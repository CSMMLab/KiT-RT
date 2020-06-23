#include "cubic2dspline.h"

Cubic2DSpline::Cubic2DSpline( const Vector& x, const Vector& y, const Matrix& data ) : _x( x ), _y( y ) { _data = addGhostLayers( data ); }
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
    for( int i = -1; i < 3; i++ )
        for( int j = -1; j < 3; j++ ) points[i + 1][j + 1] = _data( xId + j, yId + i );

    // first interpolate in x-direction and store the results on which the final interpolation in y will be done
    double interpolationBuffer[4];
    interpolationBuffer[0] = interpolate1D( points[0], x );
    interpolationBuffer[1] = interpolate1D( points[1], x );
    interpolationBuffer[2] = interpolate1D( points[2], x );
    interpolationBuffer[3] = interpolate1D( points[3], x );
    return interpolate1D( interpolationBuffer, y );
}

Matrix Cubic2DSpline::addGhostLayers( const Matrix& m ) {
    Matrix paddedMatrix( m.rows() + 4, m.columns() + 4, 0.0 );
    auto subMat = submatrix( paddedMatrix, 2u, 2u, m.rows(), m.columns() );
    subMat      = m;

    // x-direction
    for( unsigned i = 0; i < m.rows(); ++i ) {
        paddedMatrix( i, 0 )            = m( i, 0 );
        paddedMatrix( i, 1 )            = m( i, 0 );
        paddedMatrix( i, m.rows() )     = m( i, m.rows() - 1 );
        paddedMatrix( i, m.rows() + 1 ) = m( i, m.rows() - 1 );
    }

    // y-direction
    for( unsigned i = 0; i < m.columns(); ++i ) {
        paddedMatrix( 0, i )               = m( 0, i );
        paddedMatrix( 1, i )               = m( 0, i );
        paddedMatrix( m.columns(), i )     = m( m.rows() - 1, i );
        paddedMatrix( m.columns() + 1, i ) = m( m.rows() - 1, i );
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
