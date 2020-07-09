#include "catch.hpp"
#include "quadratures/quadraturebase.h"
#include "settings/globalconstants.h"

#include <vector>

std::vector<QUAD_NAME> quadraturenames         = { QUAD_MonteCarlo, QUAD_GaussLegendreTensorized, QUAD_LevelSymmetric, QUAD_Lebedev, QUAD_LDFESA };
std::vector<std::vector<int>> quadratureorders = {
    { 4, 5, 6, 7 },                            // Monte Carlo
    { 4, 6, 8, 10 },                           // Gauss Legendre
    { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 },    // Available Orders for LevelSymmetric
    { 3,  5,  7,  9,  11, 13, 15, 17, 19, 21, 23,  25,  27,  29,  31,  35,
      41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131 },    // Available orders for Lebedev
    { 1, 2, 3 }                                                                  // Available Orders for LDFESA
};

bool approxequal( double a, double b, bool lowAccuracy = false ) {
    double tol = 1e-12;              // For computed quadrature weights
    if( lowAccuracy ) tol = 5e-5;    // Mainly for Lookup Quadratures
    return std::abs( a - b ) < tol;
}

TEST_CASE( "Quadrature weights sum to 4*pi.", "[quadrature]" ) {
    bool lowAccuracyTesting = false;
    for( auto quadraturename : quadraturenames ) {
        lowAccuracyTesting = false;
        if( quadraturename == QUAD_GaussLegendreTensorized || quadraturename == QUAD_LevelSymmetric || quadraturename == QUAD_Lebedev ||
            quadraturename == QUAD_LDFESA )
            lowAccuracyTesting = true;

        for( auto quadratureorder : quadratureorders[quadraturename] ) {
            QuadratureBase* Q = QuadratureBase::CreateQuadrature( quadraturename, quadratureorder );
            if( !approxequal( Q->SumUpWeights(), 4 * M_PI, lowAccuracyTesting ) ) {
                printf( "Quadrature %d at order %d . Error : %.15f  (low accuracy testing was set to %d) \n",
                        quadraturename,
                        quadratureorder,
                        std::abs( Q->SumUpWeights() - 4 * M_PI ),
                        lowAccuracyTesting );
                printf( "Computed result %.15f \n", Q->SumUpWeights() );
            }
            REQUIRE( approxequal( Q->SumUpWeights(), 4 * M_PI, lowAccuracyTesting ) );
        }
    }
}

TEST_CASE( "Quadrature points are on the unit sphere.", "[quadrature]" ) {
    bool lowAccuracyTesting = false;
    for( auto quadraturename : quadraturenames ) {
        lowAccuracyTesting = false;
        if( quadraturename == QUAD_GaussLegendreTensorized || quadraturename == QUAD_LevelSymmetric || quadraturename == QUAD_Lebedev ||
            quadraturename == QUAD_LDFESA )
            lowAccuracyTesting = true;

        for( auto quadratureorder : quadratureorders[quadraturename] ) {
            QuadratureBase* Q   = QuadratureBase::CreateQuadrature( quadraturename, quadratureorder );
            VectorVector points = Q->GetPoints();
            for( unsigned i = 0; i < Q->GetNq(); i++ ) {
                if( !approxequal( 1.0, norm( points[i] ), lowAccuracyTesting ) ) {
                    printf( "Quadrature %d at order %d . Errorous index: %d | Error : %.15f  (low accuracy testing was set to %d) \n",
                            quadraturename,
                            quadratureorder,
                            i,
                            std::abs( norm( points[i] ) - 1.0 ),
                            lowAccuracyTesting );
                    printf( "Computed result %.15f", norm( points[i] ) );
                }
                REQUIRE( approxequal( 1.0, norm( points[i] ), lowAccuracyTesting ) );
            }
        }
    }
}

TEST_CASE( "Nq is actually equal to the number of weights.", "[quadrature]" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders[quadraturename] ) {
            QuadratureBase* Q = QuadratureBase::CreateQuadrature( quadraturename, quadratureorder );
            REQUIRE( Q->GetNq() == size( Q->GetWeights() ) );
        }
    }
}

TEST_CASE( "Nq is actually equal to the number of points.", "[quadrature]" ) {
    for( auto quadraturename : quadraturenames ) {
        for( auto quadratureorder : quadratureorders[quadraturename] ) {
            QuadratureBase* Q = QuadratureBase::CreateQuadrature( quadraturename, quadratureorder );
            REQUIRE( Q->GetNq() == size( Q->GetPoints() ) );
        }
    }
}

double f( double x, double y, double z ) {
    return x * x + y * y + z * z;    // == 1
}

double g( double x, double y, double z ) { return x + y + z; }

TEST_CASE( "Integrate a constant function.", "[quadrature]" ) {
    bool lowAccuracyTesting = false;
    for( auto quadraturename : quadraturenames ) {
        lowAccuracyTesting = false;
        if( quadraturename == QUAD_GaussLegendreTensorized || quadraturename == QUAD_LevelSymmetric || quadraturename == QUAD_Lebedev ||
            quadraturename == QUAD_LDFESA )
            lowAccuracyTesting = true;

        for( auto quadratureorder : quadratureorders[quadraturename] ) {
            QuadratureBase* Q = QuadratureBase::CreateQuadrature( quadraturename, quadratureorder );
            if( !approxequal( Q->Integrate( f ), 4.0 * M_PI, lowAccuracyTesting ) ) {
                printf( "Quadrature %d at order %d :  Error : %.15f (low accuracy testing was set to %d)\n",
                        quadraturename,
                        quadratureorder,
                        std::abs( Q->Integrate( f ) - 4.0 * M_PI ),
                        lowAccuracyTesting );
                printf( "Computed result %.15f", Q->Integrate( f ) );
            }
            REQUIRE( approxequal( Q->Integrate( f ), 4.0 * M_PI, lowAccuracyTesting ) );

            // non constant function
            // if( !approxequal( Q->Integrate( g ), -0.0997858, lowAccuracyTesting ) ) {
            //     printf( "Quadrature %d at order %d :  Error : %.15f (low accuracy testing was set to %d)\n",
            //             quadraturename,
            //             quadratureorder,
            //             std::abs( Q->Integrate( g ) + 0.0997858 ),
            //             lowAccuracyTesting );
            //     printf( "Computed result %.15f", Q->Integrate( g ) );
            // }
            // REQUIRE( approxequal( Q->Integrate( g ), -0.0997858, lowAccuracyTesting ) );
        }
    }
}
