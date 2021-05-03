#include <numeric>

#include "catch.hpp"
#include "common/config.h"
#include "common/globalconstants.h"
#include "common/io.h"
#include "common/mesh.h"

TEST_CASE( "unit mesh tests", "[mesh]" ) {
    std::string config_file_name = std::string( TESTS_PATH ) + "input/unit_tests/common/unit_mesh.cfg";

    Config* config = new Config( config_file_name );
    Mesh* mesh     = LoadSU2MeshFromFile( config );

    SECTION( "sum of all cell areas is equal to total domain volume" ) {
        double domainArea   = 1.0;
        double computedArea = 0.0;
        auto cellAreas      = mesh->GetCellAreas();
        for( auto area : cellAreas ) computedArea += area;
        REQUIRE( std::fabs( computedArea - domainArea ) < std::numeric_limits<double>::epsilon() );
    }

    SECTION( "neighbor and faces are sorted equally" ) {
        auto n                 = mesh->GetNormals();
        auto neighbors         = mesh->GetNeighbours();
        double eps             = 1e-7;
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < mesh->GetNumCells(); ++i ) {
            for( unsigned j = 0; j < mesh->GetNumNodesPerCell(); ++j ) {
                unsigned pos;
                unsigned nID = neighbors[i][j];
                if( nID == mesh->GetNumCells() ) continue;
                for( unsigned k = 0; k < neighbors[nID].size(); ++k ) {
                    if( neighbors[nID][k] == i ) pos = k;
                }
                if( blaze::l2Norm( n[i][j] + n[nID][pos] ) > eps ) errorWithinBounds = false;
            }
        }
        REQUIRE( errorWithinBounds );
    }

    SECTION( "sum over all normals yields zero" ) {
        auto n                 = mesh->GetNormals();
        double eps             = 1e-7;
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < mesh->GetNumCells(); ++i ) {
            Vector sum( 2, 0.0 );
            for( unsigned j = 0; j < mesh->GetNumNodesPerCell(); ++j ) {
                sum += n[i][j];
            }
            if( blaze::l2Norm( sum ) > eps ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }

    SECTION( "mesh does not have any unassigned faces" ) {
        auto neighbors         = mesh->GetNeighbours();
        auto boundaryType      = mesh->GetBoundaryTypes();
        bool noUnassignedFaces = true;
        for( unsigned i = 0; i < mesh->GetNumCells(); ++i ) {
            if( !( neighbors[i].size() == mesh->GetNumNodesPerCell() ||
                   ( neighbors[i].size() < mesh->GetNumNodesPerCell() && boundaryType[i] != BOUNDARY_TYPE::NONE ) ) )
                noUnassignedFaces = false;
        }
        REQUIRE( noUnassignedFaces );
    }
}
