#include <numeric>

#include "catch.hpp"
#include "io.h"

TEST_CASE( "convert image to grayscale matrix and generate suitable mesh", "[image I/O]" ) {
    std::string config_file_name = "../tests/input/linesource.cfg";

    Config* config = new Config( config_file_name );    // just to setup the logger

    std::string testImage = "../tests/input/phantom.png";
    std::string testMesh  = "../result/test.su2";
    Matrix gsImage        = createSU2MeshFromImage( testImage, testMesh );
    REQUIRE( std::filesystem::exists( testMesh ) );
    REQUIRE( blaze::max( gsImage ) > 0 );
    REQUIRE( gsImage.rows() > 0 );
    REQUIRE( gsImage.columns() > 0 );
}
