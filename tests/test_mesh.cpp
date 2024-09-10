#include <numeric>

#include "catch.hpp"
#include "common/config.hpp"
#include "common/globalconstants.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include <filesystem>
#include <iostream>

TEST_CASE( "unit mesh tests", "[mesh]" ) {
    std::string config_file_name = std::string( TESTS_PATH ) + "input/unit_tests/common/unit_mesh.cfg";

    Config* config = new Config( config_file_name );
    config->SetForcedConnectivity( true );
    Mesh* mesh = LoadSU2MeshFromFile( config );

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
                if( blaze::l2Norm( n[i][j] + n[nID][pos] ) > eps ) {
                    errorWithinBounds = false;
                    std::cout << "neighbor coordinate error: " << blaze::l2Norm( n[i][j] + n[nID][pos] ) << "\n";
                }
            }
        }
        REQUIRE( errorWithinBounds );
    }

    SECTION( "neighbor and faces are indexed consistently" ) {

        auto neighbors     = mesh->GetNeighbours();
        auto cellMidPoints = mesh->GetCellMidPoints();
        auto faceMidPoints = mesh->GetInterfaceMidPoints();
        bool isConsistent  = true;

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
            if( blaze::l2Norm( sum ) > eps ) {
                errorWithinBounds = false;

                std::cout << blaze::l2Norm( sum ) << "\n";
            }
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
    /*
    SECTION( "connectivity file is consistent with su2 file" ) {
        // Connectivity
        std::string connectivityFile = config->GetMeshFile();
        size_t lastDotIndex          = connectivityFile.find_last_of( '.' );
        connectivityFile             = connectivityFile.substr( 0, lastDotIndex );
        connectivityFile += ".con";

        if( !std::filesystem::exists( connectivityFile ) ) {
            REQUIRE( false );    // File should be written by the mesh creation
        }
        else {
            config->SetForcedConnectivity( false );
            Mesh* mesh2 = LoadSU2MeshFromFile( config );

            // Check cell number
            REQUIRE( mesh2->GetNumCells() == mesh->GetNumCells() );
            REQUIRE( mesh2->GetNumNodes() == mesh->GetNumNodes() );

            // Resize the outer vector to have nCells elements
            bool neighborOK  = true;
            bool midpointsOK = true;
            bool normalsOK   = true;
            bool boundaryOK  = true;
            double eps       = 1e-10;

            for( unsigned i = 0; i < mesh2->GetNumCells(); ++i ) {
                for( unsigned j = 0; j < mesh2->GetNumNodesPerCell(); ++j ) {
                    if( mesh2->GetNeighbours()[i][j] != mesh2->GetNeighbours()[i][j] ) {
                        neighborOK = false;
                        std::cout << "neighbor ID missmatch at index " << i << " " << j << "\n";
                    }
                    if( mesh->GetInterfaceMidPoints()[i][j].size() != mesh2->GetInterfaceMidPoints()[i][j].size() ) {
                        std::cout << "Vector size missmatch at index " << i << " " << j << " with error \n";
                        std::cout << mesh->GetInterfaceMidPoints()[i][j].size() << "\n" << mesh2->GetInterfaceMidPoints()[i][j].size() << "\n";
                    }
                    if( blaze::l2Norm( mesh->GetInterfaceMidPoints()[i][j] - mesh2->GetInterfaceMidPoints()[i][j] ) > eps ) {
                        midpointsOK = false;
                        std::cout << "midpoint missmatch at index " << i << " " << j << " with error "
                                  << blaze::l2Norm( mesh->GetInterfaceMidPoints()[i][j] - mesh2->GetInterfaceMidPoints()[i][j] ) << "\n";
                    }
                    if( blaze::l2Norm( mesh->GetNormals()[i][j] - mesh2->GetNormals()[i][j] ) > eps ) {
                        normalsOK = false;
                        std::cout << "normal missmatch at index " << i << " " << j << " with error "
                                  << blaze::l2Norm( mesh->GetInterfaceMidPoints()[i][j] - mesh2->GetInterfaceMidPoints()[i][j] ) << "\n";
                    }
                }
                if( mesh2->GetBoundaryTypes()[i] != mesh2->GetBoundaryTypes()[i] ) {
                    boundaryOK = false;
                    std::cout << "Boundary ID missmatch at index " << i << "\n";
                }
            }
            REQUIRE( neighborOK );
            REQUIRE( midpointsOK );
            REQUIRE( normalsOK );
            REQUIRE( boundaryOK );
        }
    }
    */
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
