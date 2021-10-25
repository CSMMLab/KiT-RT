#include "catch.hpp"
#include "toolboxes/interpolation.hpp"

#include <iostream>

Vector testFunction1( const Vector& x ) {
    Vector res( x.size() );
    for( unsigned i = 0; i < x.size(); ++i ) {
        res[i] = x[i];
    }
    return res;
}

Vector testFunction2( const Vector& x ) {
    Vector res( x.size() );
    for( unsigned i = 0; i < x.size(); ++i ) {
        res[i] = std::exp( x[i] );
    }
    return res;
}

Vector testFunction3( const Vector& x ) {
    Vector res( x.size() );
    for( unsigned i = 0; i < x.size(); ++i ) {
        res[i] = x[i] * x[i];
    }
    return res;
}

TEST_CASE( "interpolation tests", "[interpolation]" ) {
    SECTION( "1D - linear" ) {
        Vector x   = blaze::linspace( 10, -2.0, 2.0 );
        Vector y   = testFunction1( x );
        Vector xq  = blaze::linspace( 15, -1.0, 1.0 );
        Vector ref = testFunction1( xq );
        Interpolation interp( x, y, Interpolation::linear );
        Vector res = interp( xq );
        REQUIRE( res.size() == ref.size() );
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < res.size(); ++i ) {
            if( std::fabs( res[i] - ref[i] ) > 1e-6 ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }

    SECTION( "1D - loglinear" ) {
        Vector x   = blaze::linspace( 10, 0.5, 3.0 );
        Vector y   = testFunction2( x );
        Vector xq  = blaze::linspace( 15, 1.0, 2.5 );
        Vector ref = testFunction2( xq );
        Interpolation interp( x, y, Interpolation::loglinear );
        Vector res = interp( xq );
        REQUIRE( res.size() == ref.size() );
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < res.size(); ++i ) {
            if( std::fabs( res[i] - ref[i] ) > 1e-6 ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }

    SECTION( "1D - cubic" ) {
        Vector x   = blaze::linspace( 30, -2.0, 2.0 );
        Vector y   = testFunction3( x );
        Vector xq  = blaze::linspace( 10, -1.0, 1.0 );
        Vector ref = testFunction3( xq );
        Interpolation interp( x, y, Interpolation::cubic );
        Vector res = interp( xq );
        REQUIRE( res.size() == ref.size() );
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < res.size(); ++i ) {
            if( std::fabs( res[i] - ref[i] ) > 1e-6 ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }

    SECTION( "2D - cubic " ) {
        Vector x = blaze::linspace( 30, -2.0, 2.0 );
        Vector y = blaze::linspace( 30, -2.0, 2.0 );
        Matrix data( 30, 30, 0.0 );
        for( unsigned i = 0; i < x.size(); ++i ) {
            for( unsigned j = 0; j < y.size(); ++j ) {
                data( i, j ) = std::exp( -( x[i] / 2 + y[j] / 2 ) );
            }
        }
        Interpolation interp( x, y, data );
        Vector xq = blaze::linspace( 10, -1.5, 0.5 );
        Vector yq = blaze::linspace( 10, -1.5, 0.5 );
        Matrix ref( 10, 10, 0.0 );
        Matrix res( 10, 10, 0.0 );
        for( unsigned i = 0; i < xq.size(); ++i ) {
            for( unsigned j = 0; j < yq.size(); ++j ) {
                ref( i, j ) = std::exp( -( xq[i] / 2 + yq[j] / 2 ) );
                res( i, j ) = interp( xq[i], yq[j] );
            }
        }
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < xq.size(); ++i ) {
            for( unsigned j = 0; j < yq.size(); ++j ) {
                if( std::fabs( res( i, j ) - ref( i, j ) ) > 1e-3 ) errorWithinBounds = false;    // TODO: check low error tolerance
            }
        }
        REQUIRE( errorWithinBounds );
    }
}
