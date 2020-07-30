#include <mpi.h>

#include "io.h"
#include "solvers/solverbase.h"

#include "settings/config.h"

int main( int argc, char** argv ) {
    MPI_Init( &argc, &argv );

    std::string filename = ParseArguments( argc, argv );

    // CD  Load Settings from File
    Config* config = new Config( filename );

    // Print input file and run info to file
    PrintLogHeader( filename );

    // Build solver
    Solver* solver = Solver::Create( config );

    // Run solver and export
    solver->Solve();
    solver->Save();

    MPI_Finalize();
    return EXIT_SUCCESS;
}
