#include "catch.hpp"

#include "quadrature.h"

bool approxequal( double a, double b ) {
    double tol = 1e-15;
    return abs( a - b ) < tol;
}

TEST_CASE( "Quadrature weights sum to 4*pi.", "WHAT TO PUT HERE?" ) {
    for( int order = 4; order < 8; order++ ) {
        Quadrature* Q = Quadrature::CreateQuadrature( "montecarlo", order );
        REQUIRE( approxequal( 4 * M_PI, Q->SumUpWeights() ) );
    }
}

TEST_CASE( "Quadrature points are on the unit sphere.", "WHAT TO PUT HERE?" ) {}
