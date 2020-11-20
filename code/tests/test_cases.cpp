#include <cstdio>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

#include "catch.hpp"

#include "common/config.h"
#include "solvers/solverbase.h"

using vtkUnstructuredGridReaderSP = vtkSmartPointer<vtkUnstructuredGridReader>;

std::vector<double> readVTKFile( std::string filename ) {
    auto reader = vtkUnstructuredGridReaderSP::New();
    reader->SetFileName( filename.c_str() );
    reader->ReadAllScalarsOn();
    reader->ReadAllVectorsOn();
    reader->Update();

    auto grid         = reader->GetOutput();
    unsigned numCells = grid->GetNumberOfCells();
    auto cellData     = grid->GetCellData()->GetArray( 0 );

    std::vector<double> data( numCells, 0.0 );

    for( unsigned i = 0; i < numCells; ++i ) {
        data[i] = cellData->GetTuple1( static_cast<int>( i ) );
    }

    return data;
}

// --- Validation Tests Solver ---
TEST_CASE( "checkerboard_SN", "[validation_tests]" ) {
    std::string config_file_name = std::string( TESTS_PATH ) + "input/checkerboard_SN.cfg";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->PrintVolumeOutput();

    auto test      = readVTKFile( std::string( TESTS_PATH ) + "../result/rtsn_test_checkerboard_SN.vtk" );
    auto reference = readVTKFile( std::string( TESTS_PATH ) + "input/checkerboard_SN_reference.vtk" );

    double eps = 1e-3;
    REQUIRE( test.size() == reference.size() );
    for( unsigned i = 0; i < test.size(); ++i ) {
        REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
    }
}

TEST_CASE( "checkerboard_PN", "[validation_tests]" ) {
    std::string config_file_name = std::string( TESTS_PATH ) + "input/checkerboard_PN.cfg";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->PrintVolumeOutput();

    auto test      = readVTKFile( std::string( TESTS_PATH ) + "../result/rtsn_test_checkerboard_PN.vtk" );
    auto reference = readVTKFile( std::string( TESTS_PATH ) + "input/checkerboard_PN_reference.vtk" );

    double eps = 1e-3;
    REQUIRE( test.size() == reference.size() );
    for( unsigned i = 0; i < test.size(); ++i ) {
        REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
    }
}

TEST_CASE( "checkerboard_MN", "[validation_tests]" ) {
    std::string config_file_name = std::string( TESTS_PATH ) + "input/checkerboard_MN.cfg";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->PrintVolumeOutput();

    auto test      = readVTKFile( std::string( TESTS_PATH ) + "../result/rtsn_test_checkerboard_MN.vtk" );
    auto reference = readVTKFile( std::string( TESTS_PATH ) + "input/checkerboard_MN_reference.vtk" );

    double eps = 1e-3;
    REQUIRE( test.size() == reference.size() );
    for( unsigned i = 0; i < test.size(); ++i ) {
        REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
    }
}

TEST_CASE( "linesource_SN", "[validation_tests]" ) {
    std::string config_file_name = std::string( TESTS_PATH ) + "input/linesource_SN.cfg";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->PrintVolumeOutput();

    auto test      = readVTKFile( std::string( TESTS_PATH ) + "../result/rtsn_test_linesource_SN.vtk" );
    auto reference = readVTKFile( std::string( TESTS_PATH ) + "input/linesource_SN_reference.vtk" );

    double eps = 1e-3;
    REQUIRE( test.size() == reference.size() );
    for( unsigned i = 0; i < test.size(); ++i ) {
        REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
    }
}

TEST_CASE( "linesource_PN", "[validation_tests]" ) {
    std::string config_file_name = std::string( TESTS_PATH ) + "input/linesource_PN.cfg";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->PrintVolumeOutput();

    auto test      = readVTKFile( std::string( TESTS_PATH ) + "../result/rtsn_test_linesource_PN.vtk" );
    auto reference = readVTKFile( std::string( TESTS_PATH ) + "input/linesource_PN_reference.vtk" );

    double eps = 1e-3;

    REQUIRE( test.size() == reference.size() );
    for( unsigned i = 0; i < test.size(); ++i ) {
        REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
    }
}

TEST_CASE( "linesource_MN", "[validation_tests]" ) {

    SECTION( "Quadratic Entropy" ) {
        std::string config_file_name = std::string( TESTS_PATH ) + "input/linesource_MN_Quad.cfg";

        Config* config = new Config( config_file_name );
        Solver* solver = Solver::Create( config );
        solver->Solve();
        solver->PrintVolumeOutput();

        auto test      = readVTKFile( std::string( TESTS_PATH ) + "../result/rtsn_test_linesource_MN_Quad.vtk" );
        auto reference = readVTKFile( std::string( TESTS_PATH ) + "input/linesource_MN_Quad_reference.vtk" );

        double eps = 1e-3;

        REQUIRE( test.size() == reference.size() );
        for( unsigned i = 0; i < test.size(); ++i ) {
            REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
        }
    }

    SECTION( "Maxwell Boltzmann Entropy" ) {
        std::string config_file_name = std::string( TESTS_PATH ) + "input/linesource_MN_MB.cfg";

        Config* config = new Config( config_file_name );
        Solver* solver = Solver::Create( config );
        solver->Solve();
        solver->PrintVolumeOutput();

        auto test      = readVTKFile( std::string( TESTS_PATH ) + "../result/rtsn_test_linesource_MN_MB.vtk" );
        auto reference = readVTKFile( std::string( TESTS_PATH ) + "input/linesource_MN_MB_reference.vtk" );

        double eps = 1e-3;
        REQUIRE( test.size() == reference.size() );
        for( unsigned i = 0; i < test.size(); ++i ) {
            REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
        }
    }
}

// --- Validation Tests Output ---
TEST_CASE( "screen_output", "[output]" ) {

    spdlog::drop_all();    // Make sure to write in own logging file

    std::string config_file_name       = std::string( TESTS_PATH ) + "input/validate_logger.cfg";
    std::string screenLoggerReference  = std::string( TESTS_PATH ) + "input/validate_screen_output_reference";
    std::string screenLogger           = std::string( TESTS_PATH ) + "../result/logs/validate_screen_output";
    std::string historyLoggerReference = std::string( TESTS_PATH ) + "input/validate_screen_output_csv_reference";
    std::string historyLogger          = std::string( TESTS_PATH ) + "../result/logs/validate_screen_output_csv";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    // --- Read and validate logger ---
    std::ifstream screenLoggerReferenceStream( screenLoggerReference );
    std::ifstream screenLoggerStream( screenLogger );
    std::ifstream historyLoggerReferenceStream( historyLoggerReference );
    std::ifstream historyLoggerStream( historyLogger );

    std::string line, lineRef;
    bool lineValid;
    while( !screenLoggerReferenceStream.eof() && !screenLoggerStream.eof() ) {
        std::getline( screenLoggerReferenceStream, lineRef );
        std::getline( screenLoggerStream, line );

        lineValid = lineRef.compare( line ) == 0;
        if( !lineValid ) {
            std::cout << lineRef << "\n" << line << "\n";
        }
        REQUIRE( lineValid );
    }
    bool eqLen = screenLoggerReferenceStream.eof() && screenLoggerStream.eof();
    if( !eqLen ) {
        std::cout << "Files of unequal length!\n";
    }
    REQUIRE( eqLen );    // Files must be of same length

    while( !historyLoggerReferenceStream.eof() && !historyLoggerStream.eof() ) {
        std::getline( historyLoggerReferenceStream, lineRef );
        std::getline( historyLoggerStream, line );
        lineValid = lineRef.compare( line ) == 0;
        REQUIRE( lineValid );
    }
    eqLen = historyLoggerReferenceStream.eof() && historyLoggerStream.eof();
    if( !eqLen ) {
        std::cout << "Files of unequal length!\n";
    }
    REQUIRE( eqLen );    // Files must be of same length
}
