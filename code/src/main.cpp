#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

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

    MPI_Finalize();

    //    // Initialize the Python Interpreter
    //    std::string pyPath = RTSN_PYTHON_PATH;
    //
    //    if( !Py_IsInitialized() ) {
    //        Py_InitializeEx( 0 );
    //        if( !Py_IsInitialized() ) {
    //            ErrorMessages::Error( "Python init failed!", CURRENT_FUNCTION );
    //        }
    //        PyRun_SimpleString( ( "import sys\nsys.path.append('" + pyPath + "')" ).c_str() );
    //    }
    //
    //    PyObject *pArgs, *pReturn, *pModule, *pFunc;
    //    PyArrayObject* np_ret;
    //
    //    std::string moduleName = "callNN";
    //    pModule                = PyImport_ImportModule( moduleName.c_str() );
    //    if( !pModule ) {
    //        PyErr_Print();
    //        ErrorMessages::Error( "'" + moduleName + "' can not be imported!", CURRENT_FUNCTION );
    //    }
    //
    //    // pFunc = PyObject_GetAttrString( pModule, "initialize_network" );
    //    // if( !pFunc || !PyCallable_Check( pFunc ) ) {
    //    //     PyErr_Print();
    //    //     Py_DecRef( pModule );
    //    //     Py_DecRef( pFunc );
    //    //     ErrorMessages::Error( "'initialize_network' is null or not callable!", CURRENT_FUNCTION );
    //    // }
    //    //
    //    // // pArgs = PyTuple_New();                  // No arguments ==> need empty tuple as arguments
    //    // PyObject_CallObject( pFunc, NULL );    // No return value
    //
    //    pFunc = PyObject_GetAttrString( pModule, "call_network" );
    //    if( !pFunc || !PyCallable_Check( pFunc ) ) {
    //        PyErr_Print();
    //        Py_DecRef( pModule );
    //        Py_DecRef( pFunc );
    //        ErrorMessages::Error( "'call_network' is null or not callable!", CURRENT_FUNCTION );
    //    }
    //
    //    // -----------------
    //    //    std::string imageName   = "testimageName";
    //    //    std::string SU2Filename = "SU2Filenametest";
    //    //
    //    //    auto image = PyUnicode_FromString( imageName.c_str() );
    //
    //    double solutionOneCell[4] = { 0.0, 0.1, 0.2, 0.3 };
    //
    //    const long int dims[1] = { 4 };
    //
    //    import_array();
    //    PyObject* inputArray = PyArray_SimpleNewFromData( 1, dims, NPY_DOUBLE, (void*)solutionOneCell );
    //
    //    pArgs = PyTuple_New( 1 );
    //    PyTuple_SetItem( pArgs, 0, reinterpret_cast<PyObject*>( inputArray ) );
    //
    //    // -------------
    //    pReturn = PyObject_CallObject( pFunc, pArgs );    // PyObject
    //
    //    np_ret = reinterpret_cast<PyArrayObject*>( pReturn );    // Cast from PyObject to PyArrayObject
    //
    //    // size_t m{ static_cast<size_t>( PyArray_SHAPE( np_ret )[0] ) };
    //    // size_t n{ static_cast<size_t>( PyArray_SHAPE( np_ret )[1] ) };
    //    double* c_out = reinterpret_cast<double*>( PyArray_DATA( np_ret ) );
    //
    //    for( unsigned i = 0; i < 4; i++ ) {
    //        std::cout << c_out[i] << std::endl;
    //    }
    //
    //    // Finalizing
    //    Py_DecRef( pFunc );
    //    Py_DecRef( pModule );
    //    Py_DECREF( np_ret );
    //
    //    // Finish the Python Interpreter
    //
    //    Py_Finalize();
    //
    return 0;
}
