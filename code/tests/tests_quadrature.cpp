#include "catch.hpp"
#include "quadrature.h"

vector<std::string> quadraturenames = {"montecarlo"};

bool approxequal( double a, double b ) {
    double tol = 1e-15;
    return abs( a - b ) < tol;
}

TEST_CASE( "Quadrature weights sum to 4*pi.", "WHAT TO PUT HERE?" ) {
    for( auto quadraturename : quadraturenames ) {
        for( int order = 4; order < 8; order++ ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, order );
            REQUIRE( approxequal( 4 * M_PI, Q->SumUpWeights() ) );
        }
    }
}

TEST_CASE( "Quadrature points are on the unit sphere.", "WHAT TO PUT HERE?" ) {
    for( auto quadraturename : quadraturenames ) {
        for( int order = 4; order < 8; order++ ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, order );

            blaze::DynamicVector<blaze::DynamicVector<double>> points = Q->GetPoints();
            for( int i = 0; i < Q->GetNq(); i++ ) {
                REQUIRE( approxequal( 1.0, norm( points[i] ) ) );
            }
        }
    }
}

TEST_CASE( "Nq is actually equal to the weight vector length.", "WHAT TO PUT HERE?" ) {
    for( auto quadraturename : quadraturenames ) {
        for( int order = 4; order < 8; order++ ) {
            // Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, order );
            // REQUIRE( Q->GetNq() == size( Q->GetWeights() ) );
        }
    }
}
