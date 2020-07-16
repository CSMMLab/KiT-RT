#include "spline.h"
#include "toolboxes/errormessages.h"
#include <blaze/math/lapack/posv.h>

Spline::Spline() : m_left( first_deriv ), m_right( first_deriv ), m_left_value( 0.0 ), m_right_value( 0.0 ), m_force_linear_extrapolation( false ) {}

void Spline::set_boundary( Spline::bd_type left, double left_value, Spline::bd_type right, double right_value, bool force_linear_extrapolation ) {
    m_left                       = left;
    m_right                      = right;
    m_left_value                 = left_value;
    m_right_value                = right_value;
    m_force_linear_extrapolation = force_linear_extrapolation;
}

void Spline::set_points( const std::vector<double>& x, const std::vector<double>& y, bool cubic_Spline ) {
    Vector xn( x.size(), 0.0 );
    Vector yn( y.size(), 0.0 );
    for( unsigned i = 0; i < x.size(); ++i ) {
        xn[i] = x[i];
        yn[i] = y[i];
    }
    this->set_points( xn, yn, cubic_Spline );
}

void Spline::set_points( const Vector& x, const Vector& y, bool cubic_Spline ) {
    if( x.size() != y.size() ) ErrorMessages::Error( "Vectors are of unequal length!", CURRENT_FUNCTION );
    m_x   = x;
    m_y   = y;
    int n = x.size();
    for( int i = 0; i < n - 1; i++ ) {
        if( !( m_x[i] < m_x[i + 1] ) ) ErrorMessages::Error( "x is not sorted ascendingly!", CURRENT_FUNCTION );
    }

    if( cubic_Spline ) {
        Matrix A( n, n, 0.0 );    // TODO: should be a sparse matrix!
        Vector rhs( n, 0.0 );
        for( int i = 1; i < n - 1; i++ ) {
            A( i, i - 1 ) = 1.0 / 3.0 * ( x[i] - x[i - 1] );
            A( i, i )     = 2.0 / 3.0 * ( x[i + 1] - x[i - 1] );
            A( i, i + 1 ) = 1.0 / 3.0 * ( x[i + 1] - x[i] );
            rhs[i]        = ( y[i + 1] - y[i] ) / ( x[i + 1] - x[i] ) - ( y[i] - y[i - 1] ) / ( x[i] - x[i - 1] );
        }
        if( m_left == Spline::second_deriv ) {
            A( 0, 0 ) = 2.0;
            A( 0, 1 ) = 0.0;
            rhs[0]    = m_left_value;
        }
        else if( m_left == Spline::first_deriv ) {
            A( 0, 0 ) = 2.0 * ( x[1] - x[0] );
            A( 0, 1 ) = 1.0 * ( x[1] - x[0] );
            rhs[0]    = 3.0 * ( ( y[1] - y[0] ) / ( x[1] - x[0] ) - m_left_value );
        }
        else {
            ErrorMessages::Error( "Invalid bd_type", CURRENT_FUNCTION );
        }
        if( m_right == Spline::second_deriv ) {
            A( n - 1, n - 1 ) = 2.0;
            A( n - 1, n - 2 ) = 0.0;
            rhs[n - 1]        = m_right_value;
        }
        else if( m_right == Spline::first_deriv ) {
            A( n - 1, n - 1 ) = 2.0 * ( x[n - 1] - x[n - 2] );
            A( n - 1, n - 2 ) = 1.0 * ( x[n - 1] - x[n - 2] );
            rhs[n - 1]        = 3.0 * ( m_right_value - ( y[n - 1] - y[n - 2] ) / ( x[n - 1] - x[n - 2] ) );
        }
        else {
            ErrorMessages::Error( "Invalid bd_type", CURRENT_FUNCTION );
        }

        m_b = rhs;
        blaze::posv( A, m_b, 'U' );

        m_a.resize( n );
        m_c.resize( n );
        for( int i = 0; i < n - 1; i++ ) {
            m_a[i] = 1.0 / 3.0 * ( m_b[i + 1] - m_b[i] ) / ( x[i + 1] - x[i] );
            m_c[i] = ( y[i + 1] - y[i] ) / ( x[i + 1] - x[i] ) - 1.0 / 3.0 * ( 2.0 * m_b[i] + m_b[i + 1] ) * ( x[i + 1] - x[i] );
        }
    }
    else {
        m_a.resize( n );
        m_b.resize( n );
        m_c.resize( n );
        for( int i = 0; i < n - 1; i++ ) {
            m_a[i] = 0.0;
            m_b[i] = 0.0;
            m_c[i] = ( m_y[i + 1] - m_y[i] ) / ( m_x[i + 1] - m_x[i] );
        }
    }

    m_b0 = ( m_force_linear_extrapolation == false ) ? m_b[0] : 0.0;
    m_c0 = m_c[0];

    double h   = x[n - 1] - x[n - 2];
    m_a[n - 1] = 0.0;
    m_c[n - 1] = 3.0 * m_a[n - 2] * h * h + 2.0 * m_b[n - 2] * h + m_c[n - 2];
    if( m_force_linear_extrapolation == true ) m_b[n - 1] = 0.0;
}

double Spline::operator()( double x ) const {
    size_t n = m_x.size();
    Vector::ConstIterator it;
    it      = std::lower_bound( m_x.begin(), m_x.end(), x );
    int idx = std::max( int( it - m_x.begin() ) - 1, 0 );

    double h = x - m_x[idx];
    double interpol;
    if( x < m_x[0] ) {
        interpol = ( m_b0 * h + m_c0 ) * h + m_y[0];
    }
    else if( x > m_x[n - 1] ) {
        interpol = ( m_b[n - 1] * h + m_c[n - 1] ) * h + m_y[n - 1];
    }
    else {
        interpol = ( ( m_a[idx] * h + m_b[idx] ) * h + m_c[idx] ) * h + m_y[idx];
    }
    return interpol;
}
