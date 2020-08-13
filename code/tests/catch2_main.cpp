#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include <Python.h>
#include <mpi.h>

int main( int argc, char** argv ) {
    MPI_Init( &argc, &argv );

    wchar_t* program = Py_DecodeLocale( argv[0], NULL );
    Py_SetProgramName( program );
    const int result = Catch::Session().run( argc, argv );
    if( Py_IsInitialized() ) Py_Finalize();

    MPI_Finalize();
    return result;
}
