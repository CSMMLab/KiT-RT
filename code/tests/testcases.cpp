#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>

#include "catch.hpp"
#include "mesh.h"
#include "settings/config.h"
#include "solver.h"

using vtkUnstructuredGridReaderSP = vtkSmartPointer<vtkUnstructuredGridReader>;

std::vector<double> readVTKFile( std::string filename ) {
    auto reader = vtkUnstructuredGridReaderSP::New();
    reader->SetFileName( filename.c_str() );
    reader->ReadAllScalarsOn();
    reader->ReadAllVectorsOn();
    reader->Update();

    auto grid     = reader->GetOutput();
    auto cellData = grid->GetCellData()->GetArray( 0 );

    std::vector<double> data( cellData->GetSize(), 0.0 );

    for( unsigned i = 0; i < cellData->GetSize(); ++i ) {
        data[i] = cellData->GetTuple1( static_cast<int>( i ) );
    }

    return data;
}

TEST_CASE( "checkerboard", "[testcases]" ) {
    char config_file_name[MAX_STRING_SIZE] = "../tests/input/checkerboard.cfg";

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

TEST_CASE( "linesource", "[testcases]" ) {
    char config_file_name[MAX_STRING_SIZE] = "../tests/input/linesource.cfg";

    Config* config = new Config( config_file_name );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->Save();

    auto test      = readVTKFile( "../result/rtsn_test_linesource.vtk" );
    auto reference = readVTKFile( "../tests/input/linesource_reference.vtk" );

    double eps = 1e-3;
    REQUIRE( test.size() == reference.size() );
    for( unsigned i = 0; i < test.size(); ++i ) {
        REQUIRE( std::fabs( test[i] - reference[i] ) < eps );
    }
}
