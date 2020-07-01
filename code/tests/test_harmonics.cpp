
#include "catch.hpp"
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
    case_file.close();
}
