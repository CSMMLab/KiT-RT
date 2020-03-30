#include <mpi.h>

#include "io.h"
#include "quadrature.h"
#include "typedef.h"

int main( int argc, char** argv ) {

    MPI_Init( &argc, &argv );
    std::string inputFile = ParseArguments( argc, argv );
    Settings* settings    = ReadInputFile( inputFile );
    InitLogger( settings->GetLogDir(), spdlog::level::info, spdlog::level::info );
    PrintLogHeader( settings->GetInputFile() );

    // playground and demo
    Quadrature* Q = Quadrature::CreateQuadrature( "montecarlo", 10 );
    Q->PrintWeights();

    auto log = spdlog::get( "event" );
    log->info( "this is a print function to terminal and a logfile simultaneously" );
    log->error( "also has {0} different log types", 4 );

    Vector foo( 10, 1.0 );            // blaze vector (see typedef.h)
    std::cout << foo << std::endl;    // is printable

    log->info( BLAZE_CACHE_SIZE );
    log->info( BLAZE_USE_SHARED_MEMORY_PARALLELIZATION );

    return EXIT_SUCCESS;
}
