
#include "toolboxes/errormessages.h"

#include <mpi.h>

#include "io.h"
#include "solvers/solverbase.h"

#include "settings/config.h"

#include "solvers/sphericalharmonics.h"
#include <fstream>
#include <iostream>
#include <string>

#include <stdio.h>
using namespace std;
// ----

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

    // Finalize programm

    MPI_Finalize();

    return 0;
}
