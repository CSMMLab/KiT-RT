#include <mpi.h>

#include "io.h"

int main( int argc, char** argv ) {

    MPI_Init( &argc, &argv );
    std::string inputFile = ParseArguments( argc, argv );
    Settings* settings    = ReadInputFile( inputFile );
    InitLogger( settings->GetLogDir(), spdlog::level::info, spdlog::level::info );
    PrintLogHeader( settings->GetInputFile() );

#include "mesh.h"
    Mesh* m = LoadSU2MeshFromFile( settings );
    std::cout << &m << std::endl;

    return EXIT_SUCCESS;
}
