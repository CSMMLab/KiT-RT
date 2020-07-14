#include <Python.h>
#include <mpi.h>

#include "io.h"
#include "solvers/solverbase.h"

#include "settings/config.h"

#include "solvers/sphericalharmonics.h"
#include <fstream>
#include <iostream>
#include <string>

// ----
#include "optimizers/optimizerbase.h"
#include "quadratures/qgausslegendretensorized.h"
#include "quadratures/qmontecarlo.h"
#include "solvers/sphericalharmonics.h"

double testFunc( double my, double phi ) { return my * my + phi; }

double testFunc2( double x, double y, double z ) { return x + y + z; }

int main( int argc, char** argv ) {
    MPI_Init( &argc, &argv );
    //    wchar_t* program = Py_DecodeLocale( argv[0], NULL );
    //    Py_SetProgramName( program );

    std::string filename = ParseArguments( argc, argv );

    // CD  Load Settings from File
    Config* config = new Config( filename );

    // Test the physics reader
    // Physics testPhysic();
    // testPhysic.ReadENDL( "ENDL_H.txt" );

    // Print input file and run info to file
    PrintLogHeader( filename );

    // Build solver
    Solver* solver = Solver::Create( config );

    // Run solver and export
    solver->Solve();
    solver->Save();

    //  if( Py_IsInitialized() ) Py_Finalize();
    MPI_Finalize();
    return EXIT_SUCCESS;
}
