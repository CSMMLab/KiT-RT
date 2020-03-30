#include "catch.hpp"
#include "quadrature.h"

vector<std::string> quadraturenames = {"montecarlo"};
vector<int> quadratureorders        = {4, 5, 6, 7};

bool approxequal( double a, double b ) {
    double tol = 1e-15;
    return abs( a - b ) < tol;
}

TEST_CASE( "Quadrature weights sum to 4*pi.", "WHAT TO PUT HERE?" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            REQUIRE( approxequal( 4 * M_PI, Q->SumUpWeights() ) );
        }
    }
}

TEST_CASE( "Quadrature points are on the unit sphere.", "WHAT TO PUT HERE?" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );

            blaze::DynamicVector<blaze::DynamicVector<double>> points = Q->GetPoints();
            for( int i = 0; i < Q->GetNq(); i++ ) {
                REQUIRE( approxequal( 1.0, norm( points[i] ) ) );
            }
        }
    }
}

/*
TEST_CASE( "Nq is actually equal to the number of weights.", "WHAT TO PUT HERE?" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            REQUIRE( Q->GetNq() == size( Q->GetWeights() ) );
        }
    }
}

TEST_CASE( "Nq is actually equal to the number of points.", "WHAT TO PUT HERE?" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            REQUIRE( Q->GetNq() == size( Q->GetPoints() ) );
        }
    }
}
*/
