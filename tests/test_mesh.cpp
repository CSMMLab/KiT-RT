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

    SECTION( "neighbor and faces are indexed consistently" ) {
        auto neighbors = mesh->GetNeighbours();
        auto cellMidPoints = mesh->GetCellMidPoints();
        auto faceMidPoints = mesh->GetInterfaceMidPoints();

        bool isConsistent = true;
        for( unsigned i = 0; i < mesh->GetNumCells(); ++i ) {
            for( unsigned j = 0; j < mesh->GetNumNodesPerCell(); ++j ) {
                unsigned nID = neighbors[i][j];
                if( nID == mesh->GetNumCells() ) continue;

                auto fx = faceMidPoints[i][j][0];
                auto fy = faceMidPoints[i][j][1];

                if( fx < std::min( cellMidPoints[i][0], cellMidPoints[nID][0] ) ) isConsistent = false;
                if( fx > std::max( cellMidPoints[i][0], cellMidPoints[nID][0] ) ) isConsistent = false;
                if( fy < std::min( cellMidPoints[i][1], cellMidPoints[nID][1] ) ) isConsistent = false;
                if( fy > std::max( cellMidPoints[i][1], cellMidPoints[nID][1] ) ) isConsistent = false;
            }
        }

        REQUIRE( isConsistent );
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

TEST_CASE( "reconstruction tests", "[mesh]" ) {

    SECTION( "ensure correct Gauss theorem" ) {
        std::string config_file_name = std::string( TESTS_PATH ) + "input/unit_tests/common/unit_mesh.cfg";

        Config* config = new Config( config_file_name );
        Mesh* mesh     = LoadSU2MeshFromFile( config );

        int numCells           = mesh->GetNumCells();
        int nq                 = 5;
        auto cellMidPoints     = mesh->GetCellMidPoints();
        auto cellBoundaryTypes = mesh->GetBoundaryTypes();
        double eps             = 1e-6;

        VectorVector u( numCells, Vector( nq, 0.0 ) );

        for( int j = 0; j < numCells; ++j ) {
            u[j][0] = 1.0;    // constant
            // u[j][1] = cellMidPoints[j][0] - cellMidPoints[j][1];         // linear function x-y
            // u[j][2] = 2 * cellMidPoints[j][0] - cellMidPoints[j][1];     // linear function 2x-y
            // u[j][3] = -2 * cellMidPoints[j][0] - cellMidPoints[j][1];    // linear function -2x-y
            // u[j][4] = cellMidPoints[j][0] + cellMidPoints[j][1];         // linear function x+y
        }

        VectorVector dux( numCells, Vector( nq, 0.0 ) );
        VectorVector duy( numCells, Vector( nq, 0.0 ) );

        mesh->ComputeSlopes( nq, dux, duy, u );    // no limiter, Gauss theorem reconstruction
        bool isPass = true;

        for( int j = 0; j < numCells; ++j ) {
            if( cellBoundaryTypes[j] != 2 ) continue;
            // linear function x+y
            if( abs( dux[j][0] ) > eps || abs( duy[j][0] ) > eps ) {    // constant
                std::cout << j << " " << 0 << " : " << abs( dux[j][0] - 1.0 ) << " " << abs( duy[j][0] - 1.0 ) << std::endl;
                isPass = false;
            }
            // if( abs( dux[j][1] - 1.0 ) > eps || abs( duy[j][1] + 1.0 ) > eps ) {    // linear function x-y
            //    // std::cout << j << " " << 1 << " : " << abs( dux[j][1] - 1.0 ) << " " << abs( duy[j][1] + 1.0 ) << std::endl;
            //    isPass = false;
            //}
            // if( abs( dux[j][2] - 2.0 ) > eps || abs( duy[j][2] + 1.0 ) > eps ) {    // linear function 2x-y
            //    // std::cout << j << " " << 2 << " : " << abs( dux[j][2] - 2.0 ) << " " << abs( duy[j][2] + 1.0 ) << std::endl;
            //    isPass = false;
            //}
            // if( abs( dux[j][3] + 2.0 ) > eps || abs( duy[j][3] + 1.0 ) > eps ) {    // linear function -2x-y
            //    // std::cout << j << " " << 3 << " : " << abs( dux[j][3] + 2.0 ) << " " << abs( duy[j][3] + 1.0 ) << std::endl;
            //    isPass = false;
            //}
            // if( abs( dux[j][4] ) > eps || abs( duy[j][4] ) > eps ) {    // linear function x+y
            //    // std::cout << j << " " << 4 << " : " << ( dux[j][4] ) << " " << ( duy[j][4] ) << std::endl;
            //    isPass = false;
            //}
        }
        REQUIRE( isPass );
    }
    /*SECTION( "reconstruct correct divergence" ) {
        mesh->ReconstructSlopesU( nq, dux, duy, u );    // VK limiter
        bool isPass = true;
        for( int k = 0; k < nq; ++k ) {
            for( int j = 0; j < numCells; ++j ) {
                if( cellBoundaryTypes[j] != 2 ) continue;
                if( abs( dux[j][k] - 1.0 ) > 0.2 || abs( duy[j][k] - 1.0 ) > 0.2 ) {
                    std::cout << j << " " << dux[j][k] << " " << duy[j][k] << std::endl;
                    isPass = false;
                }
            }
        }
        REQUIRE( isPass );
    }*/
}
