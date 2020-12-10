#include "catch.hpp"
#include "common/config.h"
#include "common/globalconstants.h"
#include "quadratures/quadraturebase.h"

#include <vector>

std::vector<QUAD_NAME> quadraturenames = {
    QUAD_MonteCarlo, QUAD_GaussLegendreTensorized, QUAD_GaussLegendre1D, QUAD_LevelSymmetric, QUAD_Lebedev, QUAD_LDFESA, QUAD_Product };

std::vector<std::vector<int>> quadratureorders = {
    { 4, 5, 6, 7 },                            // Monte Carlo
    { 4, 6, 8, 10 },                           // Gauss Legendre
    { 4, 6, 8, 10 },                           // Gauss Legendre 1D
    { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 },    // Available Orders for LevelSymmetric
    { 3,  5,  7,  9,  11, 13, 15, 17, 19, 21, 23,  25,  27,  29,  31,  35,
      41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131 },    // Available orders for Lebedev
    { 1, 2, 3 },                                                                 // Available Orders for LDFESA
    { 4, 6, 8, 10 }                                                              // Available Orders for Product
};

bool approxequal( double a, double b, bool lowAccuracy = false ) {
    double tol = 1e-12;              // For computed quadrature weights
    if( lowAccuracy ) tol = 5e-5;    // Mainly for Lookup Quadratures
    return std::abs( a - b ) < tol;
}

double f( double x, double y, double z ) {
    return x * x + y * y + z * z;    // == 1
}

double sin( double x, double /*y*/, double /*z*/ ) { return sin( x ); }

TEST_CASE( "Quadrature Tests", "[quadrature]" ) {
    std::string filename = std::string( TESTS_PATH ) + "input/unit_tests/quadratures/unit_quadrature.cfg";

    // Load Settings from File
    Config* config = new Config( filename );

    SECTION( "Quadrature weights sum to 4*pi.", "[quadrature]" ) {

        bool lowAccuracyTesting = false;
        for( auto quadraturename : quadraturenames ) {
            // Set quadName
            config->SetQuadName( quadraturename );

            lowAccuracyTesting = false;
            if( quadraturename == QUAD_GaussLegendreTensorized || quadraturename == QUAD_GaussLegendre1D || quadraturename == QUAD_LevelSymmetric ||
                quadraturename == QUAD_Lebedev || quadraturename == QUAD_LDFESA || quadraturename == QUAD_Product )
                lowAccuracyTesting = true;

            for( auto quadratureorder : quadratureorders[quadraturename] ) {
                // Set quadOrder
                config->SetQuadOrder( quadratureorder );

                QuadratureBase* Q = QuadratureBase::Create( config );

                if( quadraturename == QUAD_GaussLegendre1D ) {
                    if( !approxequal( Q->SumUpWeights(), 2, lowAccuracyTesting ) ) {
                        printf( "Quadrature %d at order %d . Error : %.15f  (low accuracy testing was set to %d) ",
                                config->GetQuadName(),
                                config->GetQuadOrder(),
                                std::abs( Q->SumUpWeights() - 2 ),
                                lowAccuracyTesting );
                        printf( "Computed result %.15f \n", Q->SumUpWeights() );
                    }
                    REQUIRE( approxequal( Q->SumUpWeights(), 2, lowAccuracyTesting ) );    // 1D special case
                }
                else {
                    if( !approxequal( Q->SumUpWeights(), 4 * M_PI, lowAccuracyTesting ) ) {
                        printf( "Quadrature %d at order %d . Error : %.15f  (low accuracy testing was set to %d) ",
                                config->GetQuadName(),
                                config->GetQuadOrder(),
                                std::abs( Q->SumUpWeights() - 4 * M_PI ),
                                lowAccuracyTesting );
                        printf( "Computed result %.15f \n", Q->SumUpWeights() );
                    }
                    REQUIRE( approxequal( Q->SumUpWeights(), 4 * M_PI, lowAccuracyTesting ) );
                }
                // Special case for Gauss Legendre with half weights
                if( quadraturename == QUAD_GaussLegendreTensorized ) {
                    config->SetSNAllGaussPts( false );
                    QuadratureBase* Q = QuadratureBase::Create( config );
                    if( !approxequal( Q->SumUpWeights(), 4 * M_PI, lowAccuracyTesting ) ) {
                        printf( "Quadrature %d at order %d . Error : %.15f  (low accuracy testing was set to %d). Reduced number of quadrature "
                                "points used. \n",
                                config->GetQuadName(),
                                config->GetQuadOrder(),
                                std::abs( Q->SumUpWeights() - 4 * M_PI ),
                                lowAccuracyTesting );
                        printf( "Computed result %.15f \n", Q->SumUpWeights() );
                    }
                    REQUIRE( approxequal( Q->SumUpWeights(), 4 * M_PI, lowAccuracyTesting ) );
                    config->SetSNAllGaussPts( true );
                }
            }
        }
    }

    SECTION( "Quadrature points are on the unit sphere.", "[quadrature]" ) {

        bool lowAccuracyTesting = false;
        for( auto quadraturename : quadraturenames ) {
            // Set quadName
            config->SetQuadName( quadraturename );

            lowAccuracyTesting = false;
            if( quadraturename == QUAD_GaussLegendreTensorized || quadraturename == QUAD_GaussLegendre1D || quadraturename == QUAD_LevelSymmetric ||
                quadraturename == QUAD_Lebedev || quadraturename == QUAD_LDFESA )
                lowAccuracyTesting = true;

            if( quadraturename == QUAD_GaussLegendre1D ) continue;    // 1D test case not meaningful here

            for( auto quadratureorder : quadratureorders[quadraturename] ) {
                // Set quadOrder
                config->SetQuadOrder( quadratureorder );

                bool errorWithinBounds = true;

                QuadratureBase* Q   = QuadratureBase::Create( config );
                VectorVector points = Q->GetPoints();
                for( unsigned i = 0; i < Q->GetNq(); i++ ) {
                    if( !approxequal( 1.0, norm( points[i] ), lowAccuracyTesting ) ) {
                        printf( "Quadrature %d at order %d . Errorous index: %d | Error : %.15f  (low accuracy testing was set to %d). \n",
                                config->GetQuadName(),
                                config->GetQuadOrder(),
                                i,
                                std::abs( norm( points[i] ) - 1.0 ),
                                lowAccuracyTesting );
                        printf( "Computed result %.15f", norm( points[i] ) );
                    }
                    if( !approxequal( 1.0, norm( points[i] ), lowAccuracyTesting ) ) errorWithinBounds = false;
                }

                // Special case for Gauss Legendre with half weights
                if( quadraturename == QUAD_GaussLegendreTensorized ) {
                    config->SetSNAllGaussPts( false );
                    QuadratureBase* Q = QuadratureBase::Create( config );

                    VectorVector points = Q->GetPoints();
                    for( unsigned i = 0; i < Q->GetNq(); i++ ) {
                        if( !approxequal( 1.0, norm( points[i] ), lowAccuracyTesting ) ) {
                            printf( "Quadrature %d at order %d . Errorous index: %d | Error : %.15f  (low accuracy testing was set to %d). Reduced "
                                    "number of quadrature "
                                    "points used. \n",
                                    config->GetQuadName(),
                                    config->GetQuadOrder(),
                                    i,
                                    std::abs( norm( points[i] ) - 1.0 ),
                                    lowAccuracyTesting );
                            printf( "Computed result %.15f", norm( points[i] ) );
                        }
                        if( !approxequal( 1.0, norm( points[i] ), lowAccuracyTesting ) ) errorWithinBounds = false;
                    }

                    config->SetSNAllGaussPts( true );
                }
                REQUIRE( errorWithinBounds );
            }
        }
    }

    SECTION( "Nq is actually equal to the number of weights.", "[quadrature]" ) {

        for( auto quadraturename : quadraturenames ) {
            // Set quadName
            config->SetQuadName( quadraturename );

            for( auto quadratureorder : quadratureorders[quadraturename] ) {
                // Set quadOrder
                config->SetQuadOrder( quadratureorder );

                QuadratureBase* Q = QuadratureBase::Create( config );
                REQUIRE( Q->GetNq() == size( Q->GetWeights() ) );

                // Special case for Gauss Legendre with half weights

                if( quadraturename == QUAD_GaussLegendreTensorized ) {
                    config->SetSNAllGaussPts( false );
                    QuadratureBase* Q = QuadratureBase::Create( config );
                    REQUIRE( Q->GetNq() == size( Q->GetWeights() ) );
                    config->SetSNAllGaussPts( true );
                }
            }
        }
    }

    SECTION( "Nq is actually equal to the number of points.", "[quadrature]" ) {

        for( auto quadraturename : quadraturenames ) {
            // Set quadName
            config->SetQuadName( quadraturename );

            for( auto quadratureorder : quadratureorders[quadraturename] ) {
                // Set quadOrder
                config->SetQuadOrder( quadratureorder );

                QuadratureBase* Q = QuadratureBase::Create( config );
                REQUIRE( Q->GetNq() == size( Q->GetPoints() ) );

                // Special case for Gauss Legendre with half weights
                if( quadraturename == QUAD_GaussLegendreTensorized ) {
                    config->SetSNAllGaussPts( false );
                    QuadratureBase* Q = QuadratureBase::Create( config );
                    REQUIRE( Q->GetNq() == size( Q->GetPoints() ) );
                    config->SetSNAllGaussPts( true );
                }
            }
        }
    }

    SECTION( "Integrate a constant function.", "[quadrature]" ) {

        bool lowAccuracyTesting = false;
        for( auto quadraturename : quadraturenames ) {
            // Set quadName
            config->SetQuadName( quadraturename );

            lowAccuracyTesting = false;
            if( quadraturename == QUAD_GaussLegendreTensorized || quadraturename == QUAD_GaussLegendre1D || quadraturename == QUAD_LevelSymmetric ||
                quadraturename == QUAD_Lebedev || quadraturename == QUAD_LDFESA )
                lowAccuracyTesting = true;

            for( auto quadratureorder : quadratureorders[quadraturename] ) {
                // Set quadOrder
                config->SetQuadOrder( quadratureorder );

                QuadratureBase* Q = QuadratureBase::Create( config );

                if( quadraturename == QUAD_GaussLegendre1D ) {
                    if( !approxequal( Q->Integrate( sin ), 0, lowAccuracyTesting ) ) {
                        printf( "Quadrature %d at order %d :  Error : %.15f (low accuracy testing was set to %d)\n",
                                config->GetQuadName(),
                                config->GetQuadOrder(),
                                std::abs( Q->Integrate( sin ) ),
                                lowAccuracyTesting );
                        printf( "Computed result %.15f", Q->Integrate( sin ) );
                    }
                    REQUIRE( approxequal( Q->Integrate( sin ), 0, lowAccuracyTesting ) );    // 1D special case
                }
                else {
                    if( !approxequal( Q->Integrate( f ), 4.0 * M_PI, lowAccuracyTesting ) ) {
                        printf( "Quadrature %d at order %d :  Error : %.15f (low accuracy testing was set to %d)\n",
                                config->GetQuadName(),
                                config->GetQuadOrder(),
                                std::abs( Q->Integrate( f ) - 4.0 * M_PI ),
                                lowAccuracyTesting );
                        printf( "Computed result %.15f", Q->Integrate( f ) );
                    }
                    REQUIRE( approxequal( Q->Integrate( f ), 4.0 * M_PI, lowAccuracyTesting ) );
                }

                // Special case for Gauss Legendre with half weights
                if( quadraturename == QUAD_GaussLegendreTensorized ) {
                    config->SetSNAllGaussPts( false );
                    QuadratureBase* Q = QuadratureBase::Create( config );

                    if( !approxequal( Q->Integrate( f ), 4.0 * M_PI, lowAccuracyTesting ) ) {
                        printf(
                            "Quadrature %d at order %d :  Error : %.15f (low accuracy testing was set to %d). Reduced number of quad points used\n",
                            config->GetQuadName(),
                            config->GetQuadOrder(),
                            std::abs( Q->Integrate( f ) - 4.0 * M_PI ),
                            lowAccuracyTesting );
                        printf( "Computed result %.15f", Q->Integrate( f ) );
                    }
                    REQUIRE( approxequal( Q->Integrate( f ), 4.0 * M_PI, lowAccuracyTesting ) );

                    config->SetSNAllGaussPts( true );
                }
            }
        }
    }
}
