#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL KITRT_MLOPT_ARRAY_API
#include <numpy/arrayobject.h>

#include "common/config.h"
#include "optimizers/mloptimizer.h"
#include "toolboxes/errormessages.h"

MLOptimizer::MLOptimizer( Config* settings ) : OptimizerBase( settings ) {

    initialize_python();

    // initialize network
    std::string moduleName = "callNeuralClosure";

    _pModule = PyImport_ImportModule( moduleName.c_str() );
    if( !_pModule ) {
        PyErr_Print();
        Py_DecRef( _pModule );
        ErrorMessages::Error( "'" + moduleName + "' can not be imported!", CURRENT_FUNCTION );
    }
}

MLOptimizer::~MLOptimizer() { finalize_python(); }

void MLOptimizer::Solve( Vector& lambda, Vector& u, const VectorVector& /*moments*/, unsigned /*idx_cell*/ ) {

    // Convert Vector to array
    const unsigned input_size = u.size();
    double* nn_input          = new double[u.size()];

    for( unsigned idx_sys = 0; idx_sys < input_size; idx_sys++ ) {
        nn_input[idx_sys] = u[idx_sys];
        // std::cout << nn_input[idx_sys] << ", ";
    }

    //  initialize_python();
    double* nn_output = callNetwork( input_size, nn_input );    //  nn_input;

    // std::cout << "Solution found in cell: " << idx_cell << "/8441 \n";

    for( unsigned i = 0; i < input_size; i++ ) {
        // std::cout << nn_output[i] << ", ";
        lambda[i] = nn_output[i];
    }
    //  std::cout << std::endl;
    delete[] nn_input;
}

void MLOptimizer::SolveMultiCell( VectorVector& lambda, VectorVector& u, const VectorVector& /*moments*/ ) {

    const unsigned batch_size = u.size();       // batch size = number of cells
    const unsigned sol_dim    = u[0].size();    // dimension of input vector = nTotalEntries

    const unsigned n_size = batch_size * sol_dim;    // length of input array

    // Covert input to array
    double* nn_input = new double[n_size];

    unsigned idx_input = 0;
    for( unsigned idx_cell = 0; idx_cell < batch_size; idx_cell++ ) {
        for( unsigned idx_sys = 0; idx_sys < sol_dim; idx_sys++ ) {
            nn_input[idx_input] = u[idx_cell][idx_sys];
            idx_input++;
        }
    }

    double* nn_output = callNetworkMultiCell( batch_size, sol_dim, nn_input );

    unsigned idx_output = 0;
    for( unsigned idx_cell = 0; idx_cell < batch_size; idx_cell++ ) {
        for( unsigned idx_sys = 0; idx_sys < sol_dim; idx_sys++ ) {
            lambda[idx_cell][idx_sys] = nn_output[idx_output];
            idx_output++;
        }
    }

    delete[] nn_output;
}

void MLOptimizer::init_numpy() {
    _import_array();    // Check, if this gives a mem Leak!
}

void MLOptimizer::initialize_python() {
    // Initialize the Python Interpreter
    std::string pyPath = KITRT_PYTHON_PATH;
    pyPath             = pyPath + "/../ext/neuralEntropy/python";
    if( !Py_IsInitialized() ) {

        Py_InitializeEx( 0 );
        if( !Py_IsInitialized() ) {
            ErrorMessages::Error( "Python init failed!", CURRENT_FUNCTION );
        }
        PyRun_SimpleString( ( "import sys\nsys.path.append('" + pyPath + "')" ).c_str() );
    }

    // std::cout << "Python working directory is: " << pyPath << " \n";
    init_numpy();
}

void MLOptimizer::finalize_python() {
    Py_DecRef( _pModule );
    Py_Finalize();
}

double* MLOptimizer::callNetwork( const unsigned input_size, double* nn_input ) {

    PyObject *pArgs, *pReturn, *pFunc;    // *pModule,
    PyArrayObject* np_ret;

    pFunc = PyObject_GetAttrString( _pModule, "call_network" );
    if( !pFunc || !PyCallable_Check( pFunc ) ) {
        PyErr_Print();
        Py_DecRef( _pModule );
        Py_DecRef( pFunc );
        ErrorMessages::Error( "'call_network' is null or not callable!", CURRENT_FUNCTION );
    }

    long int dims[1] = { input_size };    // Why was this const?

    PyObject* inputArray = PyArray_SimpleNewFromData( 1, dims, NPY_DOUBLE, (void*)nn_input );

    pArgs = PyTuple_New( 1 );
    PyTuple_SetItem( pArgs, 0, reinterpret_cast<PyObject*>( inputArray ) );

    // Call Python function
    pReturn = PyObject_CallObject( pFunc, pArgs );    // PyObject

    np_ret = reinterpret_cast<PyArrayObject*>( pReturn );    // Cast from PyObject to PyArrayObject

    double* nn_output = reinterpret_cast<double*>( PyArray_DATA( np_ret ) );    // Get Output

    // Finalizing
    Py_DecRef( pFunc );
    Py_DECREF( np_ret );

    return nn_output;
}

double* MLOptimizer::callNetworkMultiCell( const unsigned batch_size, const unsigned input_dim, double* nn_input ) {

    PyObject *pArgs, *pReturn, *pFunc;
    PyArrayObject* np_ret;

    pFunc = PyObject_GetAttrString( _pModule, "call_networkBatchwise" );
    if( !pFunc || !PyCallable_Check( pFunc ) ) {
        PyErr_Print();
        Py_DecRef( _pModule );
        Py_DecRef( pFunc );
        ErrorMessages::Error( "'call_network' is null or not callable!", CURRENT_FUNCTION );
    }

    long int dims[2] = { batch_size, input_dim };    // Why was this const?

    PyObject* inputArray = PyArray_SimpleNewFromData( 2, dims, NPY_DOUBLE, (void*)nn_input );

    pArgs = PyTuple_New( 1 );
    PyTuple_SetItem( pArgs, 0, reinterpret_cast<PyObject*>( inputArray ) );

    // Call Python function
    pReturn = PyObject_CallObject( pFunc, pArgs );    // PyObject

    np_ret = reinterpret_cast<PyArrayObject*>( pReturn );    // Cast from PyObject to PyArrayObject

    double* nn_output = reinterpret_cast<double*>( PyArray_DATA( np_ret ) );    // Get Output

    // Finalizing
    Py_DecRef( pFunc );
    Py_DECREF( np_ret );

    return nn_output;
}
