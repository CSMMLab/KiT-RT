#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "optimizers/mloptimizer.h"
#include "settings/config.h"
#include "toolboxes/errormessages.h"

MLOptimizer::MLOptimizer( Config* settings ) : OptimizerBase( settings ) {
    // test
    // Test
    initialize_python();
    // initialize_Network();
}

MLOptimizer::~MLOptimizer() { finalize_python(); }

void MLOptimizer::Solve( Vector& lambda, Vector& u, VectorVector& moments, unsigned idx_cell ) {

    // Convert Vector to array
    const unsigned input_size = u.size();
    double nn_input[input_size];
    double* nn_output = new double[input_size];

    for( unsigned idx_sys = 0; idx_sys < u.size(); idx_sys++ ) {
        nn_input[idx_sys] = u[idx_sys];
        // std::cout << nn_input[idx_sys] << ", ";
    }

    //  initialize_python();
    nn_output = callNetwork( input_size, nn_input );

    for( unsigned i = 0; i < input_size; i++ ) {
        // std::cout << nn_output[i] << ", ";
        lambda[i] = nn_output[i];
    }
    // std::cout << std::endl;

    // delete[] nn_output; CHECK IF THIS PRODUCES A MEM LEAK!
}

void MLOptimizer::init_numpy() {
    _import_array();    // Check, if this gives a mem Leak!
}

void MLOptimizer::initialize_python() {
    // Initialize the Python Interpreter
    std::string pyPath = RTSN_PYTHON_PATH;

    if( !Py_IsInitialized() ) {
        Py_InitializeEx( 0 );
        if( !Py_IsInitialized() ) {
            ErrorMessages::Error( "Python init failed!", CURRENT_FUNCTION );
        }
        PyRun_SimpleString( ( "import sys\nsys.path.append('" + pyPath + "')" ).c_str() );
    }
    init_numpy();
}

void MLOptimizer::finalize_python() { Py_Finalize(); }

/*
void MLOptimizer::initialize_Network() {
    PyObject *pArgs, *pModule, *pFunc;

    std::string moduleName = "callNN";
    pModule                = PyImport_ImportModule( moduleName.c_str() );
    if( !pModule ) {
        PyErr_Print();
        ErrorMessages::Error( "'" + moduleName + "' can not be imported!", CURRENT_FUNCTION );
    }

    pFunc = PyObject_GetAttrString( pModule, "initialize_network" );
    if( !pFunc || !PyCallable_Check( pFunc ) ) {
        PyErr_Print();
        Py_DecRef( pModule );
        Py_DecRef( pFunc );
        ErrorMessages::Error( "'initialize_network' is null or not callable!", CURRENT_FUNCTION );
    }

    pArgs = PyTuple_New( 0 );               // No arguments, so empty tuple
    PyObject_CallObject( pFunc, pArgs );    // PyObject
} */

double* MLOptimizer::callNetwork( const unsigned input_size, double* nn_input ) {
    PyObject *pArgs, *pReturn, *pModule, *pFunc;
    PyArrayObject* np_ret;

    std::string moduleName = "callNN";
    pModule                = PyImport_ImportModule( moduleName.c_str() );
    if( !pModule ) {
        PyErr_Print();
        ErrorMessages::Error( "'" + moduleName + "' can not be imported!", CURRENT_FUNCTION );
    }

    pFunc = PyObject_GetAttrString( pModule, "call_network" );
    if( !pFunc || !PyCallable_Check( pFunc ) ) {
        PyErr_Print();
        Py_DecRef( pModule );
        Py_DecRef( pFunc );
        ErrorMessages::Error( "'call_network' is null or not callable!", CURRENT_FUNCTION );
    }

    const long int dims[1] = { input_size };

    PyObject* inputArray = PyArray_SimpleNewFromData( 1, dims, NPY_DOUBLE, (void*)nn_input );

    pArgs = PyTuple_New( 1 );
    PyTuple_SetItem( pArgs, 0, reinterpret_cast<PyObject*>( inputArray ) );

    // Call Python function
    pReturn = PyObject_CallObject( pFunc, pArgs );    // PyObject

    np_ret = reinterpret_cast<PyArrayObject*>( pReturn );    // Cast from PyObject to PyArrayObject

    double* c_out = reinterpret_cast<double*>( PyArray_DATA( np_ret ) );    // Get Output

    // Finalizing
    Py_DecRef( pFunc );
    Py_DecRef( pModule );
    Py_DECREF( np_ret );

    return c_out;
}
