#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

#include "catch.hpp"

#include "common/config.h"
#include "solvers/solverbase.h"
#include "toolboxes/datagenerator.h"

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

TEST_CASE( "SN_SOLVER", "[validation_tests]" ) {
    std::string sn_fileDir = "input/validation_tests/SN_solver/";
    SECTION( "checkerboard" ) {
        std::string config_file_name = std::string( TESTS_PATH ) + sn_fileDir + "checkerboard_SN.cfg";

        Config* config = new Config( config_file_name );
        Solver* solver = Solver::Create( config );

        solver->Solve();
        solver->PrintVolumeOutput();

        auto test      = readVTKFile( std::string( TESTS_PATH ) + "result/rtsn_test_checkerboard_SN.vtk" );
        auto reference = readVTKFile( std::string( TESTS_PATH ) + sn_fileDir + "checkerboard_SN_reference.vtk" );

        double eps = 1e-3;
        REQUIRE( test.size() == reference.size() );
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < test.size(); ++i ) {
            if( std::fabs( test[i] - reference[i] ) > eps ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }

    SECTION( "linesource" ) {
        std::string config_file_name = std::string( TESTS_PATH ) + sn_fileDir + "linesource_SN.cfg";

        Config* config = new Config( config_file_name );
        Solver* solver = Solver::Create( config );
        solver->Solve();
        solver->PrintVolumeOutput();

        auto test      = readVTKFile( std::string( TESTS_PATH ) + "result/rtsn_test_linesource_SN.vtk" );
        auto reference = readVTKFile( std::string( TESTS_PATH ) + sn_fileDir + "linesource_SN_reference.vtk" );

        double eps = 1e-3;
        REQUIRE( test.size() == reference.size() );
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < test.size(); ++i ) {
            if( std::fabs( test[i] - reference[i] ) > eps ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }
}

TEST_CASE( "PN_SOLVER", "[validation_tests]" ) {
    std::string pn_fileDir = "input/validation_tests/PN_solver/";
    SECTION( "checkerboard" ) {
        std::string config_file_name = std::string( TESTS_PATH ) + pn_fileDir + "checkerboard_PN.cfg";

        Config* config = new Config( config_file_name );
        Solver* solver = Solver::Create( config );
        solver->Solve();
        solver->PrintVolumeOutput();

        auto test      = readVTKFile( std::string( TESTS_PATH ) + "result/rtsn_test_checkerboard_PN.vtk" );
        auto reference = readVTKFile( std::string( TESTS_PATH ) + pn_fileDir + "checkerboard_PN_reference.vtk" );

        double eps             = 1e-3;
        bool errorWithinBounds = true;
        REQUIRE( test.size() == reference.size() );
        for( unsigned i = 0; i < test.size(); ++i ) {
            if( std::fabs( test[i] - reference[i] ) > eps ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }

    SECTION( "linesource" ) {
        std::string config_file_name = std::string( TESTS_PATH ) + pn_fileDir + "linesource_PN.cfg";

        Config* config = new Config( config_file_name );
        Solver* solver = Solver::Create( config );
        solver->Solve();
        solver->PrintVolumeOutput();

        auto test      = readVTKFile( std::string( TESTS_PATH ) + "result/rtsn_test_linesource_PN.vtk" );
        auto reference = readVTKFile( std::string( TESTS_PATH ) + pn_fileDir + "linesource_PN_reference.vtk" );

        double eps             = 1e-3;
        bool errorWithinBounds = true;

        REQUIRE( test.size() == reference.size() );
        for( unsigned i = 0; i < test.size(); ++i ) {
            if( std::fabs( test[i] - reference[i] ) > eps ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }
}

TEST_CASE( "MN_SOLVER", "[validation_tests]" ) {
    std::string mn_fileDir = "input/validation_tests/MN_solver/";

    SECTION( "checkerboard" ) {
        std::string config_file_name = std::string( TESTS_PATH ) + mn_fileDir + "checkerboard_MN.cfg";

        Config* config = new Config( config_file_name );
        Solver* solver = Solver::Create( config );
        solver->Solve();
        solver->PrintVolumeOutput();

        auto test      = readVTKFile( std::string( TESTS_PATH ) + "result/rtsn_test_checkerboard_MN.vtk" );
        auto reference = readVTKFile( std::string( TESTS_PATH ) + mn_fileDir + "checkerboard_MN_reference.vtk" );

        double eps             = 1e-3;
        bool errorWithinBounds = true;

        REQUIRE( test.size() == reference.size() );
        for( unsigned i = 0; i < test.size(); ++i ) {
            if( std::fabs( test[i] - reference[i] ) > eps ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }

    SECTION( "linesource" ) {

        // --- Quadratic Entropy
        {
            std::string config_file_name = std::string( TESTS_PATH ) + mn_fileDir + "linesource_MN_Quad.cfg";

            Config* config = new Config( config_file_name );
            Solver* solver = Solver::Create( config );
            solver->Solve();
            solver->PrintVolumeOutput();

            auto test      = readVTKFile( std::string( TESTS_PATH ) + "result/rtsn_test_linesource_MN_Quad.vtk" );
            auto reference = readVTKFile( std::string( TESTS_PATH ) + mn_fileDir + "linesource_MN_Quad_reference.vtk" );

            double eps             = 1e-3;
            bool errorWithinBounds = true;

            REQUIRE( test.size() == reference.size() );
            for( unsigned i = 0; i < test.size(); ++i ) {
                if( std::fabs( test[i] - reference[i] ) > eps ) errorWithinBounds = false;
            }
            REQUIRE( errorWithinBounds );
        }

        {    // ---  Maxwell Boltzmann Entropy ---
            std::string config_file_name = std::string( TESTS_PATH ) + mn_fileDir + "linesource_MN_MB.cfg";

            Config* config = new Config( config_file_name );
            Solver* solver = Solver::Create( config );
            solver->Solve();
            solver->PrintVolumeOutput();

            auto test      = readVTKFile( std::string( TESTS_PATH ) + "result/rtsn_test_linesource_MN_MB.vtk" );
            auto reference = readVTKFile( std::string( TESTS_PATH ) + mn_fileDir + "linesource_MN_MB_reference.vtk" );

            double eps             = 1e-3;
            bool errorWithinBounds = true;

            REQUIRE( test.size() == reference.size() );
            for( unsigned i = 0; i < test.size(); ++i ) {
                if( std::fabs( test[i] - reference[i] ) > eps ) errorWithinBounds = false;
            }
            REQUIRE( errorWithinBounds );
        }
    }
}

TEST_CASE( "CSD_SN_FP_SOLVER", "[validation_tests]" ) {
    std::string csd_sn_fileDir = "input/validation_tests/CSD_SN_FP_solver/";
    SECTION( "waterphantom 1D" ) {

        std::string config_file_name = std::string( TESTS_PATH ) + csd_sn_fileDir + "waterphantom_1D.cfg";

        Config* config = new Config( config_file_name );
        Solver* solver = Solver::Create( config );
        solver->Solve();
        solver->PrintVolumeOutput();
        auto test      = readVTKFile( std::string( TESTS_PATH ) + "result/rtsn_test_waterphantom_1D_CSD_FP.vtk" );
        auto reference = readVTKFile( std::string( TESTS_PATH ) + csd_sn_fileDir + "waterphantom_1D_CSD_SN_FP_reference.vtk" );

        double eps = 1e-3;
        REQUIRE( test.size() == reference.size() );
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < test.size(); ++i ) {
            if( std::fabs( test[i] - reference[i] ) > eps ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }
}
/*
TEST_CASE( "CSD_SN_FP_2D_SOLVER", "[validation_tests]" ) {
    std::string csd_sn_fileDir = "input/validation_tests/CSD_SN_FP_2D_solver/";
    SECTION( "waterphantom 2D" ) {

        std::string config_file_name = std::string( TESTS_PATH ) + csd_sn_fileDir + "waterphantom_2D.cfg";

        Config* config = new Config( config_file_name );
        Solver* solver = Solver::Create( config );
        solver->Solve();
        solver->PrintVolumeOutput();
        auto test      = readVTKFile( std::string( TESTS_PATH ) + "result/rtsn_test_waterphantom_2D_CSD_FP.vtk" );
        auto reference = readVTKFile( std::string( TESTS_PATH ) + csd_sn_fileDir + "waterphantom_2D_CSD_SN_FP_reference.vtk" );

        double eps = 1e-3;
        REQUIRE( test.size() == reference.size() );
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < test.size(); ++i ) {
            if( std::fabs( test[i] - reference[i] ) > eps ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }
}

TEST_CASE( "CSD_SN_FP_SH_2D_SOLVER", "[validation_tests]" ) {
    std::string csd_sn_fileDir = "input/validation_tests/CSD_SN_FP_2D_solver/";
    SECTION( "waterphantom 2D" ) {

        std::string config_file_name = std::string( TESTS_PATH ) + csd_sn_fileDir + "waterphantom_2D_sh.cfg";

        Config* config = new Config( config_file_name );
        Solver* solver = Solver::Create( config );
        solver->Solve();
        solver->PrintVolumeOutput();
        // solver->Save();
        auto test      = readVTKFile( std::string( TESTS_PATH ) + "result/waterphantom_2D_CSD_SN_FP_SH.vtk" );
        auto reference = readVTKFile( std::string( TESTS_PATH ) + csd_sn_fileDir + "waterphantom_2D_CSD_SN_FP_SH_reference.vtk" );

        double eps = 1e-3;
        REQUIRE( test.size() == reference.size() );
        bool errorWithinBounds = true;
        for( unsigned i = 0; i < test.size(); ++i ) {
            if( std::fabs( test[i] - reference[i] ) > eps ) errorWithinBounds = false;
        }
        REQUIRE( errorWithinBounds );
    }
}
*/

// --- Validation Tests Output ---
void tokenize( std::string const& str, const char delim, std::vector<std::string>& out ) {
    size_t start;
    size_t end = 0;

    while( ( start = str.find_first_not_of( delim, end ) ) != std::string::npos ) {
        end = str.find( delim, start );
        out.push_back( str.substr( start, end - start ) );
    }
}

TEST_CASE( "screen_output", "[output]" ) {
    std::string out_fileDir = "input/validation_tests/output/";
    spdlog::drop_all();    // Make sure to write in own logging file

    std::string config_file_name       = std::string( TESTS_PATH ) + out_fileDir + "validate_logger.cfg";
    std::string screenLoggerReference  = std::string( TESTS_PATH ) + out_fileDir + "validate_logger_reference";
    std::string screenLogger           = std::string( TESTS_PATH ) + "result/logs/validate_logger_output";
    std::string historyLoggerReference = std::string( TESTS_PATH ) + out_fileDir + "validate_logger_csv_reference";
    std::string historyLogger          = std::string( TESTS_PATH ) + "result/logs/validate_logger_output_csv";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();

    // Force Logger to flush
    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );
    log->flush();
    logCSV->flush();
    // --- Read and validate logger ---
    std::ifstream screenLoggerReferenceStream( screenLoggerReference );
    std::ifstream screenLoggerStream( screenLogger );
    std::ifstream historyLoggerReferenceStream( historyLoggerReference );
    std::ifstream historyLoggerStream( historyLogger );

    std::string line, lineRef;
    bool lineValid;
    const char delimScreen = '|';
    const char delimHist   = ',';

    // --- Screen Logger
    while( !screenLoggerReferenceStream.eof() && !screenLoggerStream.eof() ) {
        std::getline( screenLoggerReferenceStream, lineRef );
        std::getline( screenLoggerStream, line );

        // ignore date (before first "|")
        std::vector<std::string> out;
        std::vector<std::string> outRef;
        tokenize( line, delimScreen, out );
        tokenize( lineRef, delimScreen, outRef );

        if( out.size() != outRef.size() ) std::cout << lineRef << "\n" << line << "\n";

        REQUIRE( out.size() == outRef.size() );    // Sanity check
        // token is an colum of the csv
        for( unsigned idx_token = 1; idx_token < out.size(); idx_token++ ) {    // Skip date  ==> start from 1
            lineValid = outRef[idx_token].compare( out[idx_token] ) == 0;
            if( !lineValid ) {
                std::cout << lineRef << "\n" << line << "\n";
            }
            REQUIRE( lineValid );
        }
    }
    bool eqLen = screenLoggerReferenceStream.eof() && screenLoggerStream.eof();
    if( !eqLen ) {
        std::cout << "Files of unequal length!\n";
    }
    REQUIRE( eqLen );    // Files must be of same length

    // --- History Logger
    while( !historyLoggerReferenceStream.eof() && !historyLoggerStream.eof() ) {
        std::getline( historyLoggerReferenceStream, lineRef );
        std::getline( historyLoggerStream, line );
        // ignore date (before first "|")
        std::vector<std::string> out;
        std::vector<std::string> outRef;
        tokenize( line, delimHist, out );
        tokenize( lineRef, delimHist, outRef );

        if( out.size() != outRef.size() ) std::cout << lineRef << "\n" << line << "\n";
        REQUIRE( out.size() == outRef.size() );    // sanity check

        for( unsigned idx_token = 1; idx_token < out.size(); idx_token++ ) {    // Skip date  ==> start from 1
            lineValid = outRef[idx_token].compare( out[idx_token] ) == 0;
            if( !lineValid ) {
                std::cout << lineRef << "\n" << line << "\n";
            }
            REQUIRE( lineValid );
        }
    }
    eqLen = historyLoggerReferenceStream.eof() && historyLoggerStream.eof();
    if( !eqLen ) {
        std::cout << "Files of unequal length!\n";
    }
    REQUIRE( eqLen );    // Files must be of same length
}

TEST_CASE( "Test the Data Generator", "[dataGen]" ) {
    std::string out_fileDir = "input/validation_tests/dataGenerator/";
    spdlog::drop_all();    // Make sure to write in own logging file

    std::string config_file_name       = std::string( TESTS_PATH ) + out_fileDir + "validate_dataGen.cfg";
    std::string historyLoggerReference = std::string( TESTS_PATH ) + out_fileDir + "validate_dataGen_csv_reference";
    std::string historyLogger          = std::string( TESTS_PATH ) + "result/logs/validate_dataGen_csv";

    // --- Generate Data ---
    // Load Settings from File
    Config* config = new Config( config_file_name );
    // Build Data generator
    nnDataGenerator* datagen = new nnDataGenerator( config );
    // Generate Data and export
    datagen->computeTrainingData();

    // --- Force Logger to flush
    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );
    log->flush();
    logCSV->flush();

    // --- Read and validate logger ---
    std::ifstream historyLoggerReferenceStream( historyLoggerReference );
    std::ifstream historyLoggerStream( historyLogger );

    std::string line, lineRef;
    bool lineValid;
    const char delimHist = ',';

    bool testPassed = true;
    // --- History Logger
    unsigned count = 0;
    while( !historyLoggerReferenceStream.eof() && !historyLoggerStream.eof() && count < 3 ) {
        std::getline( historyLoggerReferenceStream, lineRef );
        std::getline( historyLoggerStream, line );
        // ignore date (before first "|")
        std::vector<std::string> out;
        std::vector<std::string> outRef;
        tokenize( line, delimHist, out );
        tokenize( lineRef, delimHist, outRef );

        if( out.size() != outRef.size() ) {
            std::cout << lineRef << "\n" << line << "\n";
            testPassed = false;
            break;
        }

        for( unsigned idx_token = 1; idx_token < out.size(); idx_token++ ) {    // Skip date  ==> start from 1
            lineValid = outRef[idx_token].compare( out[idx_token] ) == 0;
            if( !lineValid ) {
                std::cout << lineRef << "\n" << line << "\n";
                testPassed = false;
                break;
            }
        }
        count++;
    }
    REQUIRE( testPassed );
}
