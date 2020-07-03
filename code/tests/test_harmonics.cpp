
#include "catch.hpp"
#include "quadratures/qgausslegendretensorized.h"
#include "solvers/sphericalharmonics.h"

#include <fstream>
#include <sstream>

TEST_CASE( "test the spherical harmonics basis computation", "[spherical harmonics]" ) {

    std::string text_line;
    std::ifstream case_file;

    SphericalHarmonics testBase( 2 );

    double my  = 0.0;
    double phi = 0.0;
    Vector values( 9, 0.0 );
    Vector result( 4, 0.0 );

    case_file.open( "harmonicBasis_reference.csv", std::ios::in );

    // Remove title line
    getline( case_file, text_line );

    SECTION( "test to reference solution" ) {

        while( getline( case_file, text_line ) ) {

            // give line to stringstream
            std::stringstream ss( text_line );

            // Read values
            ss >> my >> phi >> values[0] >> values[1] >> values[2] >> values[3] >> values[4] >> values[5] >> values[6] >> values[7] >> values[8];

            result = testBase.ComputeSphericalBasis( my, phi );

            for( unsigned idx = 0; idx < 9; idx++ ) {
                REQUIRE( std::fabs( result[idx] - values[idx] ) < std::numeric_limits<double>::epsilon() );
            }
        }
    }

    case_file.close();

    SECTION( "test orthogonality" ) {

        QGaussLegendreTensorized quad( 6 );

        double x, y, z, w;
        Vector moment = testBase.ComputeSphericalBasis( 0, 1, 0 );
        // 9 basis moments if degree = 2

        Vector results( moment.size(), 0.0 );

        for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
            x      = quad.GetPoints()[idx_quad][0];
            y      = quad.GetPoints()[idx_quad][1];
            z      = quad.GetPoints()[idx_quad][2];
            w      = quad.GetWeights()[idx_quad];
            moment = testBase.ComputeSphericalBasis( x, y, z );

            for( unsigned idx_sys = 1; idx_sys < 9; idx_sys++ ) {
                results[idx_sys] += w * moment[idx_sys - 1] * moment[idx_sys];
            }
        }
        for( unsigned idx_sys = 0; idx_sys < 9; idx_sys++ ) {
            REQUIRE( std::fabs( results[idx_sys] ) < std::numeric_limits<double>::epsilon() );
        }
    }

    SECTION( "test normality" ) {

        QGaussLegendreTensorized quad( 6 );

        double x, y, z, w;
        Vector moment = testBase.ComputeSphericalBasis( 0, 1, 0 );
        // 9 basis moments if degree = 2

        Vector results( moment.size(), 0.0 );

        for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
            x      = quad.GetPoints()[idx_quad][0];
            y      = quad.GetPoints()[idx_quad][1];
            z      = quad.GetPoints()[idx_quad][2];
            w      = quad.GetWeights()[idx_quad];
            moment = testBase.ComputeSphericalBasis( x, y, z );

            for( unsigned idx_sys = 0; idx_sys < 9; idx_sys++ ) {
                results[idx_sys] += w * moment[idx_sys] * moment[idx_sys];
            }
        }
        for( unsigned idx_sys = 0; idx_sys < 9; idx_sys++ ) {
            REQUIRE( std::fabs( results[idx_sys] - 1 ) < std::numeric_limits<double>::epsilon() );
        }
    }
}
