#include "catch.hpp"
#include "quadrature.h"

#include <vector>

std::vector<QUAD_NAME> quadraturenames = {QUAD_MonteCarlo};
std::vector<int> quadratureorders        = {4, 5, 6, 7};

bool approxequal( double a, double b ) {
    double tol = 1e-15;
    return abs( a - b ) < tol;
}

TEST_CASE( "Quadrature weights sum to 4*pi.", "[correctweightsum]" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            REQUIRE( approxequal( 4 * M_PI, Q->SumUpWeights() ) );
        }
    }
}

TEST_CASE( "Quadrature points are on the unit sphere.", "[pointsonsphere]" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders ) {
            Quadrature* Q       = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            VectorVector points = Q->GetPoints();
            for( unsigned i = 0; i < Q->GetNq(); i++ ) {
                REQUIRE( approxequal( 1.0, norm( points[i] ) ) );
            }
        }
    }
}

TEST_CASE( "Nq is actually equal to the number of weights.", "[nqequallengthweights]" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            REQUIRE( Q->GetNq() == size( Q->GetWeights() ) );
        }
    }
}

TEST_CASE( "Nq is actually equal to the number of points.", "[nqequallengthpoints]" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            REQUIRE( Q->GetNq() == size( Q->GetPoints() ) );
        }
    }
}

double f( double x, double y, double z ) {
    return x * x + y * y + z * z;    // == 1
}

TEST_CASE( "Integrate a constant function.", "[integrateconstantfunction" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            REQUIRE( approxequal( Q->Integrate( f ), 4.0 * M_PI ) );
        }
    }
}
