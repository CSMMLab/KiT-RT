#ifndef MLOPTIMIZER_H
#define MLOPTIMIZER_H

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "optimizerbase.h"

// Forward declaration

class MLOptimizer : public OptimizerBase
{
  public:
    MLOptimizer( Config* settings );

    inline ~MLOptimizer();

    void Solve( Vector& lambda, Vector& u, const VectorVector& moments, unsigned idx_cell = 0 ) override;
    void SolveMultiCell( VectorVector& lambda, VectorVector& u, const VectorVector& moments ) override;

  private:
    /*! @brief Calls the tensorflow neural network for the entropy closure
     *  @param inputDim : dimension of moment vector for a single cell
     *  @param double* nnInput : moment vector as double array
     *  @return Lagrange multiplyer alpha with size input_size
     */
    double* callNetwork( const unsigned inputDim, double* nnInput );

    /*! @brief Calls the tensorflow neural network for the entropy closure for the whole mesh
     *  @param batchSize : number of cells in the mesh ==> batchsize for the network
     *  @param inputDim : dimension of moment vector for a single cell
     *  @param double* nnInput : moment vector as double array
     *  @return Lagrange multiplyer alpha with size input_size
     */
    double* callNetworkMultiCell( const unsigned batchSize, const unsigned inputDim, double* nnInput );

    /*! @brief Initializes the Python module. Sets Path for Python, references Python module */
    void initializePython();
    /*! @brief Initilizes numpy python module. */
    void initNumpy();
    /*! @brief Calls Python Funaction to initialize the tensorflow network. */
    void initializeNetwork();
    /*! @brief Finalizes the Python module. Dereferences Python */
    void finalizePython();

    // Python members
    PyObject* _pModule;
};

#endif    // MLOPTIMIZER_H
