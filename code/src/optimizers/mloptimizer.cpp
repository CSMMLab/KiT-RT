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
}

void MLOptimizer::Solve( Vector& lambda, Vector& u, VectorVector& moments, unsigned idx_cell ) {

    // Initialize the Python Interpreter
    std ::string pyPath = RTSN_PYTHON_PATH;

    if( !Py_IsInitialized() ) {
        Py_InitializeEx( 0 );
        if( !Py_IsInitialized() ) {
            ErrorMessages::Error( "Python init failed!", CURRENT_FUNCTION );
        }
        PyRun_SimpleString( ( "import sys\nsys.path.append('" + pyPath + "')" ).c_str() );
    }

    import_array();    // Dont know what it does, but otherwise we get a segfault DOUBLE CHECK
    if( PyErr_Occurred() ) std::cout << "error!!!\n";

    // Declare Objects
    PyObject *pArgs, *pReturn, *pModule, *pFunc;
    PyArrayObject* np_ret;

    // Import module
    std::string moduleName = "callNN";
    pModule                = PyImport_ImportModule( moduleName.c_str() );
    if( !pModule ) {
        PyErr_Print();
        ErrorMessages::Error( "'" + moduleName + "' can not be imported!", CURRENT_FUNCTION );
    }

    // Call the Neural Network as a function
    pFunc = PyObject_GetAttrString( pModule, "call_network" );
    if( !pFunc || !PyCallable_Check( pFunc ) ) {
        PyErr_Print();
        Py_DecRef( pModule );
        Py_DecRef( pFunc );
        ErrorMessages::Error( "'call_network' is null or not callable!", CURRENT_FUNCTION );
    }

    // REPLACE THAT BY A PROPER CAST
    double cell_solution[4];
    const long int dims[1] = { 4 };

    for( unsigned idx_sys = 0; idx_sys < u.size(); idx_sys++ ) {
        cell_solution[idx_sys] = u[idx_sys];
        std::cout << "," << u[idx_sys] << " ";
    }
    std::cout << std::endl;

    // double solutionOneCell[4] = { 0.0, 0.1, 0.2, 0.3 };
    //
    // const long int dims[1] = { 4 };
    std::cout << "check\n";

    std::cout << "check1\n";

    // Cast solution array to a Python object
    PyObject* inputArray = PyArray_SimpleNewFromData( 1, dims, NPY_DOUBLE, (void*)cell_solution );
    std::cout << "check2\n";

    // Pack arguments for function call
    pArgs = PyTuple_New( 1 );
    PyTuple_SetItem( pArgs, 0, reinterpret_cast<PyObject*>( inputArray ) );

    // Call the Function
    pReturn = PyObject_CallObject( pFunc, pArgs );    // PyObject

    // Prepare output ArrayObject
    np_ret = reinterpret_cast<PyArrayObject*>( pReturn );    // Cast from PyObject to PyArrayObject

    // Cast ArrayObject to c-array
    double* tempArrayOutput;
    std::cout << "check\n";
    tempArrayOutput = reinterpret_cast<double*>( PyArray_DATA( np_ret ) );
    std::cout << "check\n";

    // Write Back to Vector ==> REPLACE BY PROPER CAST
    // for( unsigned idx_sys = 0; idx_sys < u.size(); idx_sys++ ) {
    //    lambda[idx_sys] = tempArrayOutput[idx_sys];
    //}

    // Finalizing
    Py_DecRef( pFunc );
    Py_DecRef( pModule );
    Py_DECREF( np_ret );

    // Finish the Python Interpreter
    Py_Finalize();

    delete[] cell_solution;
    delete[] tempArrayOutput;
}
