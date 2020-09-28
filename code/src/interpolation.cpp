#include "interpolation.h"
#include "toolboxes/errormessages.h"
#include <blaze/math/lapack/posv.h>

Interpolation::Interpolation( const std::vector<double>& x, const std::vector<double>& y, TYPE type ) : _dim( 1u ), _type( type ) {
    _x = Vector( x.size(), x.data() );
    _y = Vector( y.size(), y.data() );
    Setup();
}

Interpolation::Interpolation( const Vector& x, const Vector& y, TYPE type ) : _dim( 1u ), _x( x ), _y( y ), _type( type ) { Setup(); }

Interpolation::Interpolation( const Vector& x, const Vector& y, const Matrix& data, TYPE type )
    : _dim( 2u ), _x( x ), _y( y ), _data( data ), _type( type ) {
    Setup();
}

void Interpolation::Setup() {
    if( _dim == 1u ) {
        if( _x.size() != _y.size() ) ErrorMessages::Error( "Vectors are of unequal length!", CURRENT_FUNCTION );
        if( _type == loglinear ) _y = blaze::log( _y );

        for( unsigned i = 0; i < _x.size() - 1u; i++ ) {
            if( !( _x[i] < _x[i + 1] ) ) ErrorMessages::Error( "x is not sorted ascendingly!", CURRENT_FUNCTION );
        }
    }
    else if( _dim == 2u ) {
        if( _x.size() != _data.rows() ) ErrorMessages::Error( "x and data are of unequal length!", CURRENT_FUNCTION );
        if( _y.size() != _data.columns() ) ErrorMessages::Error( "y and data are of unequal length!", CURRENT_FUNCTION );
        for( unsigned i = 0; i < _x.size() - 1u; i++ ) {
            if( !( _x[i] < _x[i + 1] ) ) ErrorMessages::Error( "x is not sorted ascendingly!", CURRENT_FUNCTION );
        }
        for( unsigned i = 0; i < _y.size() - 1u; i++ ) {
            if( !( _y[i] < _y[i + 1] ) ) ErrorMessages::Error( "y is not sorted ascendingly!", CURRENT_FUNCTION );
        }
    }
}

double Interpolation::operator()( double x ) const {
    if( _dim != 1u ) ErrorMessages::Error( "Invalid data dimension for operator(x)!", CURRENT_FUNCTION );

    if( x < _x[0] || x > _x[_x.size() - 1u] ) {
        ErrorMessages::Error( "Extrapolation is not supported!", CURRENT_FUNCTION );
    }

    Vector::ConstIterator it = std::lower_bound( _x.begin(), _x.end(), x );

    unsigned idx = static_cast<unsigned>( std::max( int( it - _x.begin() ) - 1, 0 ) );

    if( _type == linear || _type == loglinear ) {
        double res = _y[idx] + ( _y[idx + 1] - _y[idx] ) / ( _x[idx + 1] - _x[idx] ) * ( x - _x[idx] );

        if( _type == loglinear )
            return std::exp( res );
        else
            return res;
    }
    else if( _type == cubic ) {
        double param[4] = { _y[idx - 1], _y[idx], _y[idx + 1], _y[idx + 2] };
        double t        = ( x - _x[idx] ) / ( _x[idx + 1] - _x[idx] );
        return EvalCubic1DSpline( param, t );
    }
    else {
        ErrorMessages::Error( "Invalid type!", CURRENT_FUNCTION );
        return -1.0;
    }
}

double Interpolation::operator()( double x, double y ) const {
    if( _dim != 2u ) ErrorMessages::Error( "Invalid data dimension for operator(x,y)!", CURRENT_FUNCTION );
    if( _type == cubic ) {

        unsigned xId = IndexOfClosestValue( x, _x );
        unsigned yId = IndexOfClosestValue( y, _y );

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
        double buffer[4];
        buffer[0] = EvalCubic1DSpline( points[0], t );
        buffer[1] = EvalCubic1DSpline( points[1], t );
        buffer[2] = EvalCubic1DSpline( points[2], t );
        buffer[3] = EvalCubic1DSpline( points[3], t );
        return EvalCubic1DSpline( buffer, u );
    }
    else {
        ErrorMessages::Error( "Unsupported interpolation type!", CURRENT_FUNCTION );
        return -1.0;
    }
}

Vector Interpolation::operator()( Vector v ) const {
    Vector res( v.size() );
    for( unsigned i = 0; i < v.size(); ++i ) {
        res[i] = this->operator()( v[i] );
    }
    return res;
}

std::vector<double> Interpolation::operator()( std::vector<double> v ) const {
    std::vector<double> res( v.size() );
    for( unsigned i = 0; i < v.size(); ++i ) {
        res[i] = this->operator()( v[i] );
    }
    return res;
}

inline double Interpolation::EvalCubic1DSpline( double param[4], double x ) const {
    return ( param[1] + 0.5 * x *
                            ( param[2] - param[0] +
                              x * ( 2.0 * param[0] - 5.0 * param[1] + 4.0 * param[2] - param[3] +
                                    x * ( 3.0 * ( param[1] - param[2] ) + param[3] - param[0] ) ) ) );
}

unsigned Interpolation::IndexOfClosestValue( double value, const Vector& v ) const {
    auto i = std::min_element( begin( v ), end( v ), [=]( double x, double y ) { return abs( x - value ) < abs( y - value ); } );
    return std::distance( begin( v ), i );
}
