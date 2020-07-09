#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>

#include "catch.hpp"
#include "mesh.h"
#include "settings/config.h"
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

TEST_CASE( "checkerboard_SN", "[validation tests]" ) {
    std::string config_file_name = "../tests/input/checkerboard.cfg";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->Save();

    auto test      = readVTKFile( "../result/rtsn_test_checkerboard.vtk" );
    auto reference = readVTKFile( "../tests/input/checkerboard_reference.vtk" );

    double eps = 1e-3;
    REQUIRE( test.size() == reference.size() );
    for( unsigned i = 0; i < test.size(); ++i ) {
        REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
    }
}

TEST_CASE( "linesource_SN", "[validation tests]" ) {
    std::string config_file_name = "../tests/input/linesource_SN.cfg";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->Save();

    auto test      = readVTKFile( "../result/rtsn_test_linesource_SN.vtk" );
    auto reference = readVTKFile( "../tests/input/linesource_SN_reference.vtk" );

    double eps = 1e-3;
    REQUIRE( test.size() == reference.size() );
    for( unsigned i = 0; i < test.size(); ++i ) {
        REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
    }
}

TEST_CASE( "linesource_PN", "[validation tests]" ) {
    char config_file_name[MAX_STRING_SIZE] = "../tests/input/linesource_PN.cfg";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->Save();

    auto test      = readVTKFile( "../result/rtsn_test_linesource_PN.vtk" );
    auto reference = readVTKFile( "../tests/input/linesource_PN_reference.vtk" );

    double eps = 1e-3;
    REQUIRE( test.size() == reference.size() );
    for( unsigned i = 0; i < test.size(); ++i ) {
        REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
    }
}

TEST_CASE( "linesource_MN", "[validation tests]" ) {
    char config_file_name[MAX_STRING_SIZE] = "../tests/input/linesource_MN.cfg";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->Save();

    auto test      = readVTKFile( "../result/rtsn_test_linesource_PN.vtk" );
    auto reference = readVTKFile( "../tests/input/linesource_MN_reference.vtk" );

    double eps = 1e-3;
    REQUIRE( test.size() == reference.size() );
    for( unsigned i = 0; i < test.size(); ++i ) {
        REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
    }
}
