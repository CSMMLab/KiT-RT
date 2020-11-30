#include <fstream>
#include <numeric>

#include "catch.hpp"
#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "toolboxes/interpolation.h"
#include "toolboxes/textprocessingtoolbox.h"

TEST_CASE( "convert image data to grayscale matrix", "[image I/O]" ) {
    std::string config_file_name = std::string( TESTS_PATH ) + "input/unit_tests/common/image_conversion.cfg";

    Config* config = new Config( config_file_name );    // just to init spdlog

    std::string testImage = std::string( TESTS_PATH ) + "input/unit_tests/common/mini_phantom.png";
    std::string testMesh  = config->GetMeshFile();

    Matrix gsImage = createSU2MeshFromImage( testImage, testMesh );
    gsImage.transpose();
    SECTION( "grayscale matrix" ) {
        REQUIRE( std::filesystem::exists( testMesh ) );    // mesh has been created
        REQUIRE( gsImage.rows() > 0 );                     // atleast some data is stored
        REQUIRE( gsImage.columns() > 0 );                  //
        REQUIRE( blaze::min( gsImage ) >= 0 );             // lower bound
        REQUIRE( blaze::max( gsImage ) <= 1.0 );           // upper bound

        // load reference matrix from csv file
        std::string refMatrixFile = std::string( TESTS_PATH ) + "input/unit_tests/common/phantom.csv";
        std::ifstream data( refMatrixFile );
        REQUIRE( data.is_open() );
        std::string line;
        std::vector<std::vector<double>> refMatrix;
        while( std::getline( data, line ) ) {
            std::stringstream lineStream( line );
            std::string cell;
            std::vector<double> row;
            while( std::getline( lineStream, cell, ',' ) ) {
                row.push_back( std::stod( cell ) );
            }
            refMatrix.push_back( row );
        }

        REQUIRE( refMatrix.size() == gsImage.rows() );          // equal number of rows
        REQUIRE( refMatrix[0].size() == gsImage.columns() );    // equal number of columns
        REQUIRE( std::all_of( begin( refMatrix ), end( refMatrix ), [refMatrix]( const std::vector<double>& x ) {
            return x.size() == refMatrix[0].size();
        } ) );    // consistency check if all columns of the read-in file have equal length

        bool matricesEqual = true;
        for( unsigned i = 0; i < gsImage.rows(); ++i ) {
            for( unsigned j = 0; j < gsImage.columns(); ++j ) {
                if( refMatrix[i][j] != gsImage( i, j ) ) matricesEqual = false;    // all values match
            }
        }
        REQUIRE( matricesEqual );
    }

    SECTION( "interpolation of grayscale matrix onto the generated mesh" ) {
        Mesh* mesh       = LoadSU2MeshFromFile( config );
        auto cellCenters = mesh->GetCellMidPoints();
        auto bounds      = mesh->GetBounds();

        double xMin = bounds[0].first;
        double xMax = bounds[0].second;
        double yMin = bounds[1].first;
        double yMax = bounds[1].second;

        unsigned m = gsImage.rows();
        unsigned n = gsImage.columns();

        Vector x( m ), y( n );
        for( unsigned i = 0; i < m; ++i ) x[i] = static_cast<double>( i ) / static_cast<double>( m - 1 ) * ( xMax - xMin );
        for( unsigned i = 0; i < n; ++i ) y[i] = static_cast<double>( i ) / static_cast<double>( n - 1 ) * ( yMax - yMin );

        Interpolation interp( x, y, gsImage );
        std::vector<double> result( mesh->GetNumCells(), 0.0 );
        for( unsigned i = 0; i < mesh->GetNumCells(); ++i ) {
            result[i] = std::clamp( interp( cellCenters[i][0], cellCenters[i][1] ), 0.0, 1.0 );
        }

        std::vector<std::string> fieldNames{ "CT Data" };
        std::vector<std::vector<std::string>> fieldNamesWrapper{ fieldNames };

        std::vector<std::vector<double>> scalarField( 1, result );
        std::vector<std::vector<std::vector<double>>> results{ scalarField };
        std::string outputFile = config->GetOutputFile();
        if( !TextProcessingToolbox::StringEndsWith( outputFile, ".vtk" ) ) outputFile.append( ".vtk" );
        ExportVTK( outputFile, results, fieldNamesWrapper, mesh );

        REQUIRE( std::filesystem::exists( outputFile ) );

        delete mesh;

        std::remove( outputFile.c_str() );
    }

    std::remove( testMesh.c_str() );
}
