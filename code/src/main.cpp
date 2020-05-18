#include <mpi.h>

#include "io.h"
#include "solver.h"

#include "settings/config.h"

int main( int argc, char** argv ) {

    MPI_Init( &argc, &argv );

    std::string filename = "default.cfg";
    char config_file_name[MAX_STRING_SIZE];

    filename = ParseArguments( argc, argv );

    /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
     file is specified, default.cfg is used) ---*/
    strcpy( config_file_name, filename.c_str() );

    // Load Settings from File
    Config* config = new Config( config_file_name );

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
