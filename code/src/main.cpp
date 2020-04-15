#include <mpi.h>

#include "io.h"

int main( int argc, char** argv ) {

    MPI_Init( &argc, &argv );
    std::string inputFile = ParseArguments( argc, argv );
    Settings* settings    = ReadInputFile( inputFile );
    InitLogger( settings->GetLogDir(), spdlog::level::info, spdlog::level::info );
    PrintLogHeader( settings->GetInputFile() );

#include "mesh.h"
    Mesh* m     = LoadSU2MeshFromFile( settings );
    auto colors = m->GetPartitionIDs();
    std::vector<double> dcolors( colors.size() );
    for( unsigned i = 0; i < colors.size(); ++i ) dcolors[i] = static_cast<double>( colors[i] );
    std::vector<std::vector<double>> scalarField{dcolors};
    std::vector<std::vector<std::vector<double>>> out{scalarField};
    std::vector<std::string> fieldNames{"color"};
    ExportVTK( "partitioning", out, fieldNames, settings, m );
    MPI_Finalize();
    return EXIT_SUCCESS;
}
