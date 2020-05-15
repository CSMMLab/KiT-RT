#include <mpi.h>

#include "io.h"
#include "solver.h"

#include "settings/config.h"

int main( int argc, char** argv ) {

    MPI_Init( &argc, &argv );

    /*std::string inputFile = ParseArguments( argc, argv );
    Settings* settings    = ReadInputFile( inputFile );
    InitLogger( settings->GetLogDir(), spdlog::level::info, spdlog::level::info );
    PrintLogHeader( settings->GetInputFile() );

    // build solver
    Solver* solver = Solver::Create( settings );
    solver->Solve();
    solver->Save();
    */

    std::string filename = "default.cfg";
    char config_file_name[MAX_STRING_SIZE];

    filename = ParseArguments( argc, argv );

    /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
     file is specified, default.cfg is used) ---*/
    strcpy(config_file_name, filename.c_str());

    //Load Settings from File
    Config* config = new Config(config_file_name);


    // build solver
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->Save();


    // TODO: call solver
    MPI_Finalize();
    return EXIT_SUCCESS;
}
