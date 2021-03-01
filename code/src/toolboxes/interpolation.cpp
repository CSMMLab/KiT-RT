#include "toolboxes/interpolation.h"
#include "blaze/math/dense/DenseIterator.h"        // for operator-, DenseIter...
#include "blaze/math/dense/DynamicVector.h"        // for DynamicVector<>::Con...
#include "blaze/math/expressions/DVecMapExpr.h"    // for DVecMapExpr
#include "blaze/math/expressions/Vector.h"         // for begin, end
#include "blaze/math/smp/default/DenseVector.h"    // for smpAssign
#include "toolboxes/errormessages.h"
#include <algorithm>    // for lower_bound, max
#include <iterator>     // for distance
#include <math.h>       // for log, exp
#include <stdlib.h>     // for abs

// Change Vector type to blaze for std input
Interpolation::Interpolation( const std::vector<double>& x, const std::vector<double>& y, TYPE type ) : _dim( 1u ), _type( type ) {
    _x = Vector( x.size(), x.data() );
    _y = Vector( y.size(), y.data() );
    Setup();
}

// 1D Constructor: Initialise values and call Set-up
Interpolation::Interpolation( const Vector& x, const Vector& y, TYPE type ) : _dim( 1u ), _x( x ), _y( y ), _type( type ) { Setup(); }

// 2D Constructor: Initialise values and call Set-up
Interpolation::Interpolation( const Vector& x, const Vector& y, const Matrix& data, TYPE type )
    : _dim( 2u ), _x( x ), _y( y ), _data( data ), _type( type ) {
    Setup();
}

// Set-up: check validity of input
void Interpolation::Setup() {
    // 1D
    if( _dim == 1u ) {
        // Length of tables must be equal
        if( _x.size() != _y.size() ) ErrorMessages::Error( "Vectors are of unequal length!", CURRENT_FUNCTION );
        if( _type == loglinear ) _y = blaze::log( _y );

        // Values must be sorted ascendingly (property is used later when searching tables)
        for( unsigned i = 0; i < _x.size() - 1u; i++ ) {
            if( !( _x[i] < _x[i + 1] ) ) ErrorMessages::Error( "x is not sorted ascendingly!", CURRENT_FUNCTION );
        }
    }

    // 2D
    else if( _dim == 2u ) {
        // Length of tables must be equal
        if( _x.size() != _data.rows() ) ErrorMessages::Error( "x and data are of unequal length!", CURRENT_FUNCTION );
        if( _y.size() != _data.columns() ) ErrorMessages::Error( "y and data are of unequal length!", CURRENT_FUNCTION );

        // Values must be sorted ascendingly (property is used later when searching tables)
        for( unsigned i = 0; i < _x.size() - 1u; i++ ) {
            if( !( _x[i] < _x[i + 1] ) ) ErrorMessages::Error( "x is not sorted ascendingly!", CURRENT_FUNCTION );
        }
        for( unsigned i = 0; i < _y.size() - 1u; i++ ) {
            if( !( _y[i] < _y[i + 1] ) ) ErrorMessages::Error( "y is not sorted ascendingly!", CURRENT_FUNCTION );
        }
    }
}

// 1D interpolation
double Interpolation::operator()( double x ) const {
    // Check whether 1D
    if( _dim != 1u ) ErrorMessages::Error( "Invalid data dimension for operator(x)!", CURRENT_FUNCTION );
    // x must be between min and max of table values
    if( x < _x[0] || x > _x[_x.size() - 1u] ) {
        // std::cout << x << "\t" << _x[0] << std::endl;
        // std::cout << x << "\t" << _x[_x.size() - 1u] << std::endl;
        ErrorMessages::Error( "Extrapolation is not supported!", CURRENT_FUNCTION );
    }

    // points to first value in _x that is not smaller than x (upper bound)
    Vector::ConstIterator it = std::lower_bound( _x.begin(), _x.end(), x );

    // index of the lower bound to x in _x
    unsigned idx = static_cast<unsigned>( std::max( int( it - _x.begin() ) - 1, 0 ) );

    // linear interpolation
    if( _type == linear || _type == loglinear ) {
        // interpolate between lower bound and upper bound of x in _x
        double res = _y[idx] + ( _y[idx + 1] - _y[idx] ) / ( _x[idx + 1] - _x[idx] ) * ( x - _x[idx] );

        if( _type == loglinear )
            return std::exp( res );
        else
            return res;
    }

    // cubic interpolation
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

// 2D interpolation
double Interpolation::operator()( double x, double y ) const {
    // Check whether 2D
    if( _dim != 2u ) ErrorMessages::Error( "Invalid data dimension for operator(x,y)!", CURRENT_FUNCTION );
    if( _type == cubic ) {

        // find closest values to x and y in table (lower bounds)
        int xId = IndexOfClosestValue( x, _x );
        int yId = IndexOfClosestValue( y, _y );

        // store all 16 interpolation points needed
        double points[4][4];
        for( int i = -1; i < 3; ++i ) {
            unsigned idx_y;
            idx_y = yId + i < 0 ? 0 : yId + i;
            idx_y = yId + i > static_cast<int>( _data.rows() - 1 ) ? _data.rows() - 1 : yId + i;
            for( int j = -1; j < 3; ++j ) {
                unsigned idx_x;
                idx_x = xId + j < 0 ? 0 : xId + j;
                idx_x = xId + j > static_cast<int>( _data.columns() - 1 ) ? _data.columns() - 1 : xId + j;

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

// extension of 1D interpolation to Vector input
Vector Interpolation::operator()( Vector v ) const {
    Vector res( v.size() );
    for( unsigned i = 0; i < v.size(); ++i ) {
        res[i] = this->operator()( v[i] );
    }
    return res;
}

// extension of 1D interpolation to Vector input for std::vectors
std::vector<double> Interpolation::operator()( std::vector<double> v ) const {
    std::vector<double> res( v.size() );
    for( unsigned i = 0; i < v.size(); ++i ) {
        res[i] = this->operator()( v[i] );
    }
    return res;
}

// implementation of third degree polynomial
inline double Interpolation::EvalCubic1DSpline( double param[4], double x ) const {
    return ( param[1] + 0.5 * x *
                            ( param[2] - param[0] +
                              x * ( 2.0 * param[0] - 5.0 * param[1] + 4.0 * param[2] - param[3] +
                                    x * ( 3.0 * ( param[1] - param[2] ) + param[3] - param[0] ) ) ) );
}

// find index of closes value to 'value' in vector 'v'
unsigned Interpolation::IndexOfClosestValue( double value, const Vector& v ) const {
    auto i = std::min_element( begin( v ), end( v ), [=]( double x, double y ) { return abs( x - value ) < abs( y - value ); } );
    return std::distance( begin( v ), i );
}
