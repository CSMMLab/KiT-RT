#include <mpi.h>

#include "io.h"
#include "solver.h"

int main( int argc, char** argv ) {

    MPI_Init( &argc, &argv );
    std::string inputFile = ParseArguments( argc, argv );
    Settings* settings    = ReadInputFile( inputFile );
    InitLogger( settings->GetLogDir(), spdlog::level::info, spdlog::level::info );
    PrintLogHeader( settings->GetInputFile() );

    // build solver
    Solver* solver = Solver::Create( settings );
    solver->Solve();
    solver->Save();

    // TODO: call solver
    MPI_Finalize();
    return EXIT_SUCCESS;
}
