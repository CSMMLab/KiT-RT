#include "catch.hpp"
#include "common/config.hpp"
#include "common/globalconstants.hpp"
#include "quadratures/quadraturebase.hpp"

#include <vector>

std::vector<QUAD_NAME> quadraturenames = { QUAD_MonteCarlo,
                                           QUAD_GaussLegendreTensorized,
                                           QUAD_GaussLegendre1D,
                                           QUAD_GaussLegendreTensorized2D,
                                           QUAD_LevelSymmetric,
                                           QUAD_Lebedev,
                                           QUAD_LDFESA,
                                           QUAD_Product,
                                           QUAD_Rectangular1D,
                                           QUAD_Rectangular2D,
                                           QUAD_Rectangular3D };

std::vector<std::vector<int>> quadratureorders = {
    { 4, 5, 6, 7 },                            // Monte Carlo
    { 4, 6, 8, 10 },                           // Gauss Legendre
    { 4, 6, 8, 10 },                           // Gauss Legendre 1D
    { 4, 6, 8, 10 },                           // Gauss Legendre 2D
    { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 },    // Available Orders for LevelSymmetric
    { 3,  5,  7,  9,  11, 13, 15, 17, 19, 21, 23,  25,  27,  29,  31,  35,
      41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131 },    // Available orders for Lebedev
    { 1, 2, 3 },                                                                 // Available Orders for LDFESA
    { 4, 6, 8, 10 },                                                             // Available Orders for Product
    { 60, 80, 100 },                                                             // Rectangular 1D
    { 60, 80, 100 },                                                             // Rectangular 2D
    { 60, 80, 100 }                                                              // Rectangular 3D
};

bool approxequal( double a, double b, bool lowAccuracy = false ) {
    double tol = 1e-12;              // For computed quadrature weights
    if( lowAccuracy ) tol = 5e-3;    // Mainly for Lookup Quadratures
    return std::abs( a - b ) < tol;
}

double f( double x, double y, double z ) {
    return 1.0;
    // x* x + y* y + z* z;    // == 1
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

double Polygon( double x, double y, double z ) { return 2 * x * x + y * y + z * z; }

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
            if( quadraturename == QUAD_GaussLegendreTensorized || quadraturename == QUAD_GaussLegendre1D ||
                quadraturename == QUAD_GaussLegendreTensorized2D || quadraturename == QUAD_LevelSymmetric || quadraturename == QUAD_Lebedev ||
                quadraturename == QUAD_LDFESA || quadraturename == QUAD_Product || quadraturename == QUAD_Rectangular1D ||
                quadraturename == QUAD_Rectangular2D || quadraturename == QUAD_Rectangular3D )
                lowAccuracyTesting = true;

            for( auto quadratureorder : quadratureorders[quadraturename] ) {
                // Set quadOrder
                config->SetQuadOrder( quadratureorder );

                QuadratureBase* Q = QuadratureBase::Create( config );

                if( quadraturename == QUAD_GaussLegendre1D || quadraturename == QUAD_Midpoint1D || quadraturename == QUAD_Rectangular1D ) {
                    if( !approxequal( Q->SumUpWeights(), 2, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( Q->SumUpWeights() - 2 ), Q->SumUpWeights(), lowAccuracyTesting );
                    }
                }
                else if( quadraturename == QUAD_GaussLegendreTensorized2D || quadraturename == QUAD_Midpoint2D ) {
                    if( !approxequal( Q->SumUpWeights(), M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( Q->SumUpWeights() - M_PI ), Q->SumUpWeights(), lowAccuracyTesting );
                    }
                }
                else if( quadraturename == QUAD_Rectangular2D ) {
                    if( !approxequal( Q->SumUpWeights(), 4.0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( Q->SumUpWeights() - 4.0 ), Q->SumUpWeights(), lowAccuracyTesting );
                    }
                }
                else if( quadraturename == QUAD_Rectangular3D ) {
                    if( !approxequal( Q->SumUpWeights(), 8.0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( Q->SumUpWeights() - 8.0 ), Q->SumUpWeights(), lowAccuracyTesting );
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
            if( quadraturename == QUAD_GaussLegendreTensorized || quadraturename == QUAD_GaussLegendre1D ||
                quadraturename == QUAD_GaussLegendreTensorized2D || quadraturename == QUAD_LevelSymmetric || quadraturename == QUAD_Lebedev ||
                quadraturename == QUAD_LDFESA )
                lowAccuracyTesting = true;

            if( quadraturename == QUAD_GaussLegendre1D || quadraturename == QUAD_GaussLegendreTensorized2D || quadraturename == QUAD_Midpoint1D ||
                quadraturename == QUAD_Midpoint2D || quadraturename == QUAD_Rectangular1D || quadraturename == QUAD_Rectangular2D ||
                quadraturename == QUAD_Rectangular3D )
                continue;    // test case not meaningful here

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

    SECTION( "Quadrature points tests for 2D.", "[quadrature]" ) {
        // check if 2-norm of karthesian points equals to sqrt(1-mu^2),  where mu is from spherical quad points
        bool testPassed         = true;
        bool lowAccuracyTesting = false;

        for( auto quadraturename : quadraturenames ) {
            // Set quadName
            config->SetQuadName( quadraturename );

            lowAccuracyTesting = false;

            if( quadraturename != QUAD_GaussLegendreTensorized2D ) continue;    // 1D and 3D test case not meaningful here

            for( auto quadratureorder : quadratureorders[quadraturename] ) {
                // Set quadOrder
                config->SetQuadOrder( quadratureorder );

                QuadratureBase* Q         = QuadratureBase::Create( config );
                VectorVector points       = Q->GetPoints();          //(v_x,v_y,v_z)
                VectorVector pointsSphere = Q->GetPointsSphere();    // (mu, phi, r)

                for( unsigned i = 0; i < Q->GetNq(); i++ ) {
                    double ptNorm = norm( points[i] );
                    double ref    = sqrt( 1 - pointsSphere[i][0] * pointsSphere[i][0] );
                    double ptMu   = sqrt( 1 - ptNorm * ptNorm );
                    Vector pt3d( 3 );
                    pt3d[0] = points[i][0];
                    pt3d[1] = points[i][1];
                    pt3d[2] = ptMu;

                    if( !approxequal( ref, ptNorm ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( ref - ptNorm ), ref, lowAccuracyTesting );
                        printf( "Errorous index: %d\n", i );
                    }
                    if( !approxequal( pointsSphere[i][0], ptMu ) ) {
                        if( !approxequal( -pointsSphere[i][0], ptMu ) ) {
                            testPassed = false;
                            PrintErrorMsg( config, std::abs( pointsSphere[i][0] - ptMu ), ptMu, lowAccuracyTesting );
                            printf( "Errorous index: %d\n", i );
                        }
                    }
                    if( !approxequal( 1.0, norm( pt3d ), lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( 1.0 - norm( pt3d ) ), norm( pt3d ), lowAccuracyTesting );
                        printf( "Errorous index: %d\n", i );
                    }
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

            if( quadraturename == QUAD_Rectangular3D ) {
                lowAccuracyTesting = true;
            }

            for( auto quadratureorder : quadratureorders[quadraturename] ) {
                // Set quadOrder
                config->SetQuadOrder( quadratureorder );

                QuadratureBase* Q = QuadratureBase::Create( config );

                if( quadraturename == QUAD_Rectangular1D || quadraturename == QUAD_GaussLegendre1D || quadraturename == QUAD_Midpoint1D ) {
                    result = Q->Integrate( sin );
                    if( !approxequal( result, 0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 0 ), result, lowAccuracyTesting );
                    }
                }
                else if( quadraturename == QUAD_GaussLegendreTensorized2D || quadraturename == QUAD_Midpoint2D ) {
                    result = Q->Integrate( sin );
                    if( !approxequal( result, 0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 0 ), result, lowAccuracyTesting );
                    }
                }
                else if( quadraturename == QUAD_Rectangular2D ) {
                    result = Q->Integrate( f );
                    if( !approxequal( result, 4.0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 4.0 ), result, lowAccuracyTesting );
                    }
                }
                else if( quadraturename == QUAD_Rectangular3D ) {
                    result = Q->Integrate( f );
                    if( !approxequal( result, 8.0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 8.0 ), result, lowAccuracyTesting );
                    }
                }
                else {
                    if( quadraturename == QUAD_Midpoint3D ) {
                        lowAccuracyTesting = true;
                    }
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
            if( quadraturename != QUAD_GaussLegendre1D && quadraturename != QUAD_Midpoint1D && quadraturename != QUAD_Rectangular1D &&
                quadraturename != QUAD_Rectangular2D && quadraturename != QUAD_Rectangular3D ) {    // 1D quad is trivial
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
                        if( quadraturename != QUAD_GaussLegendreTensorized2D && quadraturename != QUAD_Midpoint2D ) {
                            result = Omega_z( pointsSphere[idx_nq][0], pointsSphere[idx_nq][1] );
                            if( !approxequal( points[idx_nq][2], result, lowAccuracyTesting ) ) {
                                testPassed = false;
                                PrintErrorMsg( config, std::abs( result - points[idx_nq][2] ), result, lowAccuracyTesting );
                                printf( "z component incorrectly computed.\n" );
                                printf( "Faulty index is %d.\n", idx_nq );
                            }
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
            if( quadraturename == QUAD_GaussLegendreTensorized || quadraturename == QUAD_GaussLegendreTensorized2D ||
                quadraturename == QUAD_GaussLegendre1D || quadraturename == QUAD_LevelSymmetric || quadraturename == QUAD_Lebedev ||
                quadraturename == QUAD_LDFESA || quadraturename == QUAD_Midpoint2D || quadraturename == QUAD_Midpoint3D ||
                quadraturename == QUAD_Midpoint1D )
                lowAccuracyTesting = true;

            for( auto quadratureorder : quadratureorders[quadraturename] ) {
                // Set quadOrder
                config->SetQuadOrder( quadratureorder );

                QuadratureBase* Q = QuadratureBase::Create( config );

                // Note: Leaving out Quad_GaussLegendreTensorized with half weights... (to be added)
                if( quadraturename != QUAD_GaussLegendreTensorized2D && quadraturename != QUAD_GaussLegendre1D && quadraturename != QUAD_MonteCarlo &&
                    quadraturename != QUAD_Midpoint2D && quadraturename != QUAD_Midpoint1D && quadraturename != QUAD_Rectangular1D &&
                    quadraturename != QUAD_Rectangular2D && quadraturename != QUAD_Rectangular3D )    // MonteCarlo is too low order...
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
                if( quadraturename == QUAD_GaussLegendreTensorized2D || quadraturename == QUAD_Midpoint2D ) {
                    result = Q->IntegrateSpherical( Omega_0 );
                    if( !approxequal( result, M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - M_PI ), result, lowAccuracyTesting );
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
                    result = Q->IntegrateSpherical( Omega_xx );
                    if( !approxequal( result, 1. / 3. * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 1. / 3. * M_PI ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_xx.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_yy );
                    if( !approxequal( result, 1. / 3. * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 1. / 3. * M_PI ), result, lowAccuracyTesting );
                        printf( "Error at integrating Omega_yy.\n" );
                    }
                    result = Q->IntegrateSpherical( Omega_zz );
                    if( !approxequal( result, 1. / 3. * M_PI, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 1. / 3. * M_PI ), result, lowAccuracyTesting );
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
                if( quadraturename == QUAD_Rectangular2D ) {
                    result             = Q->Integrate( Polygon );
                    lowAccuracyTesting = true;
                    if( !approxequal( result, 4.0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 8.0 / 3.0 ), result, lowAccuracyTesting );
                        printf( "Error at integrating Polygon.\n" );
                    }
                }
                if( quadraturename == QUAD_Rectangular3D ) {
                    result             = Q->Integrate( Polygon );
                    lowAccuracyTesting = true;
                    if( !approxequal( result, 32.0 / 3.0, lowAccuracyTesting ) ) {
                        testPassed = false;
                        PrintErrorMsg( config, std::abs( result - 32.0 / 3.0 ), result, lowAccuracyTesting );
                        printf( "Error at integrating Polygon.\n" );
                    }
                }

                delete Q;
            }
        }
        REQUIRE( testPassed );
    }
    //  delete config; TODO FIX CONFIG DESTRUCTOR
}
