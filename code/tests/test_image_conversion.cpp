#include <numeric>

#include "catch.hpp"
#include "io.h"

TEST_CASE( "convert image data to grayscale matrix", "[image I/O]" ) {
    std::string config_file_name = "../tests/input/image_conversion.cfg";

    Config* config = new Config( config_file_name );    // just to init spdlog

    std::string testImage = "../tests/input/phantom.png";
    std::string testMesh  = config->GetMeshFile();

    Matrix gsImage = createSU2MeshFromImage( testImage, testMesh );
    SECTION( "grayscale matrix" ) {
        REQUIRE( std::filesystem::exists( testMesh ) );    // mesh has been created
        REQUIRE( gsImage.rows() > 0 );                     // atleast some data is stored
        REQUIRE( gsImage.columns() > 0 );                  //
        REQUIRE( blaze::min( gsImage ) >= 0 );             // lower bound
        REQUIRE( blaze::max( gsImage ) <= 1.0 );           // upper bound

        // load reference matrix from csv file
        std::string refMatrixFile = "../tests/input/phantom.csv";
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

        for( unsigned i = 0; i < gsImage.rows(); ++i ) {
            for( unsigned j = 0; j < gsImage.columns(); ++j ) {
                REQUIRE( refMatrix[i][j] == gsImage( i, j ) );    // all values match
            }
        }
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

        Vector x( m + 1 ), y( n + 1 );
        for( unsigned i = 0; i < m + 1; ++i ) {
            x[i] = static_cast<double>( i ) / static_cast<double>( m ) * ( xMax - xMin );
        }
        for( unsigned i = 0; i < n + 1; ++i ) y[i] = static_cast<double>( i ) / static_cast<double>( n ) * ( yMax - yMin );

        Cubic2DSpline interp( x, y, gsImage );
        std::vector<double> result( mesh->GetNumCells(), 0.0 );
        for( unsigned i = 0; i < mesh->GetNumCells(); ++i ) {
            result[i] = std::clamp( interp( cellCenters[i][0], cellCenters[i][1] ), 0.0, 1.0 );
        }

        std::vector<std::string> fieldNames{ "CT Data" };
        std::vector<std::vector<double>> scalarField( 1, result );
        std::vector<std::vector<std::vector<double>>> results{ scalarField };
        std::string outputFile = config->GetOutputFile();
        if( !TextProcessingToolbox::StringEndsWith( outputFile, ".vtk" ) ) outputFile.append( ".vtk" );
        ExportVTK( outputFile, results, fieldNames, mesh );

        REQUIRE( std::filesystem::exists( outputFile ) );

        delete mesh;

        // std::remove( outputFile.c_str() );
    }

    std::remove( testMesh.c_str() );
}
