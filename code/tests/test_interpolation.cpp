#include "catch.hpp"
#include "toolboxes/interpolation.h"

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
    SECTION( "linear" ) {
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

    SECTION( "loglinear" ) {
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

    SECTION( "cubic" ) {
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
}
