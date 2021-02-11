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

double Omega_0( double /*my*/, double /*phi*/ ) { return 1; }

double Omega_x( double my, double phi ) { return sqrt( 1 - my * my ) * cos( phi ); }
double Omega_y( double my, double phi ) { return sqrt( 1 - my * my ) * sin( phi ); }
double Omega_z( double my, double /*phi*/ ) { return my; }

double Omega_xx( double my, double phi ) { return Omega_x( my, phi ) * Omega_x( my, phi ); }
double Omega_yy( double my, double phi ) { return Omega_y( my, phi ) * Omega_y( my, phi ); }
double Omega_zz( double my, double phi ) { return Omega_z( my, phi ) * Omega_z( my, phi ); }
double Omega_xy( double my, double phi ) { return Omega_x( my, phi ) * Omega_y( my, phi ); }
double Omega_xz( double my, double phi ) { return Omega_x( my, phi ) * Omega_z( my, phi ); }
double Omega_yz( double my, double phi ) { return Omega_y( my, phi ) * Omega_z( my, phi ); }

void PrintErrorMsg( Config* config, double absErr, double result, bool lowAccuracyTesting ) {
    printf( "Quadrature %d at order %d . Error : %.15f  (low accuracy testing was set to %d) ",
            config->GetQuadName(),
            config->GetQuadOrder(),
            absErr,
            lowAccuracyTesting );
    printf( "Computed result %.15f \n", result );
}

TEST_CASE( "Quadrature Tests", "[quadrature]" ) {
    std::string filename = std::string( TESTS_PATH ) + "input/unit_tests/quadratures/unit_quadrature.cfg";

    // Load Settings from File
    Config* config = new Config( filename );

    SECTION( "Quadrature weights sum to 4*pi.", "[quadrature]" ) {
        bool testPassed         = true;
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
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( Q->SumUpWeights() - 2 ), Q->SumUpWeights(), lowAccuracyTesting );
                    }

                }
                else {
                    if( !approxequal( Q->SumUpWeights(), 4 * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( Q->SumUpWeights() - 4 * M_PI ), Q->SumUpWeights(), lowAccuracyTesting );
                    }

                }
                // Special case for Gauss Legendre with half weights
                if( quadraturename == QUAD_GaussLegendreTensorized ) {
                    config->SetSNAllGaussPts( false );
                    QuadratureBase* Q = QuadratureBase::Create( config );
                    if( !approxequal( Q->SumUpWeights(), 4 * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( Q->SumUpWeights() - 4 * M_PI ), Q->SumUpWeights(), lowAccuracyTesting );
                        printf( "Reduced number of quadrature was points used. \n" );
                    }

                    config->SetSNAllGaussPts( true );
                }
                delete Q;
            }
        }
        REQUIRE( testPassed );
    }

    SECTION( "Quadrature points are on the unit sphere.", "[quadrature]" ) {
        bool testPassed         = true;
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

                QuadratureBase* Q   = QuadratureBase::Create( config );
                VectorVector points = Q->GetPoints();
                for( unsigned i = 0; i < Q->GetNq(); i++ ) {
                    if( !approxequal( 1.0, norm( points[i] ), lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( norm( points[i] ) - 1.0 ), norm( points[i] ), lowAccuracyTesting );
                        printf( "Errorous index: %d\n", i );
                    }
                }

                // Special case for Gauss Legendre with half weights
                if( quadraturename == QUAD_GaussLegendreTensorized ) {
                    config->SetSNAllGaussPts( false );
                    QuadratureBase* Q = QuadratureBase::Create( config );

                    VectorVector points = Q->GetPoints();
                    for( unsigned i = 0; i < Q->GetNq(); i++ ) {
                        if( !approxequal( 1.0, norm( points[i] ), lowAccuracyTesting ) ) {
                            testPassed = false;
                            PrintErrorMsg( config, std::abs( norm( points[i] ) - 1.0 ), norm( points[i] ), lowAccuracyTesting );
                            printf( "Errorous index: %d\n", i );
                            printf( "Reduced number of quadrature was points used. \n" );
                        }
                    }

                    config->SetSNAllGaussPts( true );
                }
                delete Q;
            }
        }
        REQUIRE( testPassed );
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
                delete Q;
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
                delete Q;
            }
        }
    }

    SECTION( "Integrate a constant function.", "[quadrature]" ) {
        double result           = 0.0;
        bool testPassed         = true;
        bool lowAccuracyTesting = false;

        for( auto quadraturename : quadraturenames ) {
            // Set quadName
            config->SetQuadName( quadraturename );

            for( auto quadratureorder : quadratureorders[quadraturename] ) {
                // Set quadOrder
                config->SetQuadOrder( quadratureorder );

                QuadratureBase* Q = QuadratureBase::Create( config );

                if( quadraturename == QUAD_GaussLegendre1D ) {
                    result = Q->Integrate( sin );
                    if( !approxequal( result, 0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 0 ), result, lowAccuracyTesting );
                    }

                }
                else {
                    result = Q->Integrate( f );
                    if( !approxequal( result, 4.0 * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 4.0 * M_PI ), result, lowAccuracyTesting );
                    }

                }

                // Special case for Gauss Legendre with half weights
                if( quadraturename == QUAD_GaussLegendreTensorized ) {
                    config->SetSNAllGaussPts( false );
                    QuadratureBase* Q = QuadratureBase::Create( config );

                    result = Q->Integrate( f );
                    if( !approxequal( result, 4.0 * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 4.0 * M_PI ), result, lowAccuracyTesting );
                        printf( "Reduced number of quadrature was points used. \n" );
                    }


                    config->SetSNAllGaussPts( true );
                }
            }
        }
        REQUIRE( testPassed );
    }

    SECTION( "Consistency between carthesian and spherical coordinates", "[quadrature]" ) {
        double result           = 0.0;
        bool testPassed         = true;
        bool lowAccuracyTesting = false;

        VectorVector points;
        VectorVector pointsSphere;

        for( auto quadraturename : quadraturenames ) {
            if( quadraturename != QUAD_GaussLegendre1D ) {    // 1D quad is trivial
                config->SetQuadName( quadraturename );
                for( auto quadratureorder : quadratureorders[quadraturename] ) {
                    config->SetQuadOrder( quadratureorder );
                    QuadratureBase* Q = QuadratureBase::Create( config );
                    points            = Q->GetPoints();
                    pointsSphere      = Q->GetPointsSphere();
                    for( unsigned idx_nq = 0; idx_nq < Q->GetNq(); idx_nq++ ) {
                        result = Omega_x( pointsSphere[idx_nq][0], pointsSphere[idx_nq][1] );
                        if( !approxequal( points[idx_nq][0], result, lowAccuracyTesting ) ) {
                            testPassed = false;
                            PrintErrorMsg( config, std::abs( result - points[idx_nq][0] ), result, lowAccuracyTesting );
                            printf( "x component incorrectly computed.\n" );
                            printf( "Faulty index is %d.\n", idx_nq );
                        }
                        result = Omega_y( pointsSphere[idx_nq][0], pointsSphere[idx_nq][1] );
                        if( !approxequal( points[idx_nq][1], result, lowAccuracyTesting ) ) {
                            testPassed = false;
                            PrintErrorMsg( config, std::abs( result - points[idx_nq][1] ), result, lowAccuracyTesting );
                            printf( "y component incorrectly computed.\n" );
                            printf( "Faulty index is %d.\n", idx_nq );
                        }
                        result = Omega_z( pointsSphere[idx_nq][0], pointsSphere[idx_nq][1] );
                        if( !approxequal( points[idx_nq][2], result, lowAccuracyTesting ) ) {
                            testPassed = false;
                            PrintErrorMsg( config, std::abs( result - points[idx_nq][2] ), result, lowAccuracyTesting );
                            printf( "z component incorrectly computed.\n" );
                            printf( "Faulty index is %d.\n", idx_nq );
                        }
                    }
                    delete Q;
                }
            }
        }

        REQUIRE( testPassed );
    }

    SECTION( "Integrate polynomials", "[quadrature]" ) {

        double result           = 0.0;
        bool testPassed         = true;
        bool lowAccuracyTesting = false;

        for( auto quadraturename : quadraturenames ) {

            // Leave out Lookup Quadratures
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

                // Note: Leaving out Quad_GaussLegendreTensorized with half weights... (to be added)
                if( quadraturename != QUAD_GaussLegendre1D && quadraturename != QUAD_MonteCarlo )    // MonteCarlo is too low order...
                {
                    if( quadraturename == QUAD_LevelSymmetric && quadratureorder == 20 ) continue;    // Order 20 is somehow errorous
                    result = Q->IntegrateSpherical( Omega_0 );
                    if( !approxequal( result, 4.0 * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 4.0 * M_PI ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_0.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_x );
                    if( !approxequal( result, 0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 0 ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_x.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_y );
                    if( !approxequal( Q->IntegrateSpherical( Omega_y ), 0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 0 ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_y.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_z );
                    if( !approxequal( result, 0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 0 ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_z.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_xx );
                    if( !approxequal( result, 4. / 3. * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 4. / 3. * M_PI ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_xx.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_yy );
                    if( !approxequal( result, 4. / 3. * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 4. / 3. * M_PI ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_yy.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_zz );
                    if( !approxequal( result, 4. / 3. * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 4. / 3. * M_PI ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_zz.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_xy );
                    if( !approxequal( result, 0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 0 ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_xy.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_xz );
                    if( !approxequal( result, 0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 0 ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_xz.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_yz );
                    if( !approxequal( result, 0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 0 ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_yz.\n" );
                    }
                }
                delete Q;
            }
        }
        REQUIRE( testPassed );
    }
    //  delete config; TODO FIX CONFIG DESTRUCTOR
}
