#include <mpi.h>

#include "io.h"
#include "solver.h"

#include "settings/CConfig.h"

int main( int argc, char** argv ) {

    MPI_Init( &argc, &argv );

    char config_file_name[MAX_STRING_SIZE];

    std::string filename = ParseArguments( argc, argv );

    /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
     file is specified, default.cfg is used) ---*/
    strcpy( config_file_name, filename.c_str() );

    // Load Settings from File
    CConfig* config = new CConfig( config_file_name );

    // Init stdout and file logger
    InitLogger( config->GetLogDir(), spdlog::level::info, spdlog::level::info );

    // Print input file and run info to file
    PrintLogHeader( config_file_name );

    // Build solver
    Solver* solver = Solver::Create( config );

    // Run solver and export
    solver->Solve();
    solver->Save();

    MPI_Finalize();
    return EXIT_SUCCESS;
}
