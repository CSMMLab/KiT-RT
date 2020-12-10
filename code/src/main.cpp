/*! @file: main.cpp
 *  @brief: Main method to call the KiT-RT solver suite
 *  @author: J. Kusch, S. Schotthöfer, P. Stammer, J. Wolters, T. Xiao
 *  @version: 0.1
 */

#include <Python.h>
#include <mpi.h>
#include <string>

#include "common/config.h"
#include "common/io.h"
#include "solvers/solverbase.h"

#include "toolboxes/datagenerator.h"

int main( int argc, char** argv ) {
    MPI_Init( &argc, &argv );
    wchar_t* program = Py_DecodeLocale( argv[0], NULL );
    Py_SetProgramName( program );

    std::string filename = ParseArguments( argc, argv );

    // CD  Load Settings from File
    Config* config = new Config( filename );

    // Print input file and run info to file
    PrintLogHeader( filename );

    if( config->GetDataGeneratorMode() ) {
        // Build Data generator
        nnDataGenerator* datagen = new nnDataGenerator( config );
        // Generate Data and export
        datagen->computeTrainingData();
    }
    else {
        // Build solver
        Solver* solver = Solver::Create( config );

        // Run solver and export
        solver->Solve();
        solver->PrintVolumeOutput();
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
