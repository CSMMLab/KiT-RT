#include "interpolation.h"
#include "toolboxes/errormessages.h"
#include <blaze/math/lapack/posv.h>

Interpolation::Interpolation( const std::vector<double>& x, const std::vector<double>& y, TYPE type )
    : _left( first_deriv ), _right( first_deriv ), _type( type ), _left_value( 0.0 ), _right_value( 0.0 ), _force_linear_extrapolation( false ) {
    this->set_points( x, y, type );
}

Interpolation::Interpolation( const Vector& x, const Vector& y, TYPE type )
    : _left( first_deriv ), _right( first_deriv ), _type( type ), _left_value( 0.0 ), _right_value( 0.0 ), _force_linear_extrapolation( false ) {
    this->set_points( x, y, type );
}

void Interpolation::set_boundary( BOUNDARY left, double left_value, BOUNDARY right, double right_value, bool force_linear_extrapolation ) {
    _left                       = left;
    _right                      = right;
    _left_value                 = left_value;
    _right_value                = right_value;
    _force_linear_extrapolation = force_linear_extrapolation;
}

void Interpolation::set_points( const std::vector<double>& x, const std::vector<double>& y, TYPE type ) {
    Vector xn( x.size(), 0.0 );
    Vector yn( y.size(), 0.0 );
    for( unsigned i = 0; i < x.size(); ++i ) {
        xn[i] = x[i];
        yn[i] = y[i];
    }
    this->set_points( xn, yn, type );
}

void Interpolation::set_points( const Vector& x, const Vector& y, TYPE type ) {
    if( x.size() != y.size() ) ErrorMessages::Error( "Vectors are of unequal length!", CURRENT_FUNCTION );
    _x = x;
    if( type == loglinear )
        _y = log( y );
    else
        _y = y;
    int n = x.size();
    for( int i = 0; i < n - 1; i++ ) {
        if( !( _x[i] < _x[i + 1] ) ) ErrorMessages::Error( "x is not sorted ascendingly!", CURRENT_FUNCTION );
    }

    if( type == cubic ) {
        Matrix A( n, n, 0.0 );    // TODO: should be a sparse matrix!
        Vector rhs( n, 0.0 );
        for( int i = 1; i < n - 1; i++ ) {
            A( i, i - 1 ) = 1.0 / 3.0 * ( x[i] - x[i - 1] );
            A( i, i )     = 2.0 / 3.0 * ( x[i + 1] - x[i - 1] );
            A( i, i + 1 ) = 1.0 / 3.0 * ( x[i + 1] - x[i] );
            rhs[i]        = ( y[i + 1] - y[i] ) / ( x[i + 1] - x[i] ) - ( y[i] - y[i - 1] ) / ( x[i] - x[i - 1] );
        }
        if( _left == Interpolation::second_deriv ) {
            A( 0, 0 ) = 2.0;
            A( 0, 1 ) = 0.0;
            rhs[0]    = _left_value;
        }
        else if( _left == Interpolation::first_deriv ) {
            A( 0, 0 ) = 2.0 * ( x[1] - x[0] );
            A( 0, 1 ) = 1.0 * ( x[1] - x[0] );
            rhs[0]    = 3.0 * ( ( y[1] - y[0] ) / ( x[1] - x[0] ) - _left_value );
        }
        else {
            ErrorMessages::Error( "Invalid bd_type", CURRENT_FUNCTION );
        }
        if( _right == Interpolation::second_deriv ) {
            A( n - 1, n - 1 ) = 2.0;
            A( n - 1, n - 2 ) = 0.0;
            rhs[n - 1]        = _right_value;
        }
        else if( _right == Interpolation::first_deriv ) {
            A( n - 1, n - 1 ) = 2.0 * ( x[n - 1] - x[n - 2] );
            A( n - 1, n - 2 ) = 1.0 * ( x[n - 1] - x[n - 2] );
            rhs[n - 1]        = 3.0 * ( _right_value - ( y[n - 1] - y[n - 2] ) / ( x[n - 1] - x[n - 2] ) );
        }
        else {
            ErrorMessages::Error( "Invalid boundary type!", CURRENT_FUNCTION );
        }

        _b = rhs;
        blaze::posv( A, _b, 'U' );

        _a.resize( n );
        _c.resize( n );
        for( int i = 0; i < n - 1; i++ ) {
            _a[i] = 1.0 / 3.0 * ( _b[i + 1] - _b[i] ) / ( x[i + 1] - x[i] );
            _c[i] = ( y[i + 1] - y[i] ) / ( x[i + 1] - x[i] ) - 1.0 / 3.0 * ( 2.0 * _b[i] + _b[i + 1] ) * ( x[i + 1] - x[i] );
        }
    }
    else if( type == linear || type == loglinear ) {
        _a.resize( n );
        _b.resize( n );
        _c.resize( n );
        for( int i = 0; i < n - 1; i++ ) {
            _a[i] = 0.0;
            _b[i] = 0.0;
            _c[i] = ( _y[i + 1] - _y[i] ) / ( _x[i + 1] - _x[i] );
        }
    }
    else {
        ErrorMessages::Error( "Invalid interpolation type!", CURRENT_FUNCTION );
    }

    _b0 = ( _force_linear_extrapolation == false ) ? _b[0] : 0.0;
    _c0 = _c[0];

    double h  = x[n - 1] - x[n - 2];
    _a[n - 1] = 0.0;
    _c[n - 1] = 3.0 * _a[n - 2] * h * h + 2.0 * _b[n - 2] * h + _c[n - 2];
    if( _force_linear_extrapolation == true ) _b[n - 1] = 0.0;
}

double Interpolation::operator()( double x ) const {
    size_t n = _x.size();
    Vector::ConstIterator it;
    it      = std::lower_bound( _x.begin(), _x.end(), x );
    int idx = std::max( int( it - _x.begin() ) - 1, 0 );

    double h = x - _x[idx];
    double interpol;
    if( x < _x[0] ) {
        interpol = ( _b0 * h + _c0 ) * h + _y[0];
    }
    else if( x > _x[n - 1] ) {
        interpol = ( _b[n - 1] * h + _c[n - 1] ) * h + _y[n - 1];
    }
    else {
        interpol = ( ( _a[idx] * h + _b[idx] ) * h + _c[idx] ) * h + _y[idx];
    }
    if( _type == loglinear )
        return std::exp( interpol );
    else
        return interpol;
}

Vector Interpolation::operator()( Vector v ) const {
    Vector res( v.size() );
    for( unsigned i = 0; i < v.size(); ++i ) {
        res[i] = this->operator()( v[i] );
    }
    return res;
}
