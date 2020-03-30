#include "io.h"
#include "quadrature.h"
#include <mpi.h>

int main( int argc, char** argv ) {

    MPI_Init( &argc, &argv );
    std::string inputFile = ParseArguments( argc, argv );
    Settings* settings    = ReadInputFile( inputFile );
    InitLogger( settings->GetLogDir(), spdlog::level::info, spdlog::level::info );
    PrintLogHeader( settings->GetInputFile() );
    Quadrature* Q = Quadrature::CreateQuadrature( "montecarlo", 10 );
    Q->PrintWeights();
    return EXIT_SUCCESS;
}
