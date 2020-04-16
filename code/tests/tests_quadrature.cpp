#include "catch.hpp"
#include "../../include/quadratures/quadrature.h"
#include "option_structure.h"

#include <vector>

std::vector<QUAD_NAME> quadraturenames = {QUAD_MonteCarlo , QUAD_LevelSymmetric, QUAD_Lebedev ,QUAD_LDFESA};
std::vector<std::vector<int>> quadratureorders = {{4, 5, 6, 7}, //Monte Carlo
                                                  {},           //Gauss Legendre not working right now
                                                  {2, 4, 6 , 8, 10, 12, 14, 16, 18, 20}, //Available Orders for LevelSymmetric
                                                  {3, 5, 7 , 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53,
                                                   59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131}, //Available orders for Lebedev
                                                  {1, 2, 3}  //Available Orders for LDFESA
                                                 };

bool approxequal( double a, double b ) {
    double tol = 1e-15;//1e-5//1e-15;
    return abs( a - b ) < tol;
}

TEST_CASE( "Quadrature weights sum to 4*pi.", "[correctweightsum]" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders[quadraturename] ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            if(! approxequal( Q->SumUpWeights(), 4 * M_PI )){
                printf("Quadrature %d at order %d . Error : %.15f  \n",quadraturename, quadratureorder, abs(  Q->SumUpWeights() - 4 * M_PI  )  );
            }
            REQUIRE( approxequal( Q->SumUpWeights(), 4 * M_PI ) );
        }
    }
}

TEST_CASE( "Quadrature points are on the unit sphere.", "[pointsonsphere]" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders[quadraturename] ) {
            Quadrature* Q       = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            VectorVector points = Q->GetPoints();
            for( unsigned i = 0; i < Q->GetNq(); i++ ) {
                if(! approxequal( 1.0, norm( points[i] )  )) {
                    printf("Quadrature %d at order %d . Errorous index: %d | Error : %.15f  \n",quadraturename, quadratureorder, i,abs( norm( points[i] ) - 1.0  )  );
                }
                REQUIRE( approxequal( 1.0, norm( points[i] ) ) );
            }
        }
    }
}

TEST_CASE( "Nq is actually equal to the number of weights.", "[nqequallengthweights]" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders[quadraturename] ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            REQUIRE( Q->GetNq() == size( Q->GetWeights() ) );
        }
    }
}

TEST_CASE( "Nq is actually equal to the number of points.", "[nqequallengthpoints]" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders[quadraturename] ) {
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
        for( auto quadratureorder : quadratureorders[quadraturename] ) {
            Quadrature* Q = Quadrature::CreateQuadrature( quadraturename, quadratureorder );
            if(! approxequal( Q->Integrate( f ), 4.0 * M_PI ) ) {
                printf("Quadrature %d at order %d :  Error : %.15f \n",quadraturename, quadratureorder, abs(  Q->Integrate( f ) - 4.0 * M_PI  )  );
            }
            REQUIRE( approxequal( Q->Integrate( f ), 4.0 * M_PI ) );
        }
    }
}
