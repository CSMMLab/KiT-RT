#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include <Python.h>
#include <filesystem>
#include <mpi.h>

int main( int argc, char** argv ) {
    MPI_Init( &argc, &argv );

    wchar_t* program = Py_DecodeLocale( argv[0], NULL );
    Py_SetProgramName( program );
    const int result = Catch::Session().run( argc, argv );
    if( Py_IsInitialized() ) Py_Finalize();

    std::filesystem::remove_all( std::string( TESTS_PATH ) + "result" );

    MPI_Finalize();
    return result;
}
