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
    double* callNetwork( const unsigned input_size, double* input );
    double* callNetworkMultiCell( const unsigned batch_size, const unsigned input_dim, double* nn_input );

    void finalize_python();
    void initialize_python();
    void init_numpy();

    // Python members
    PyObject* _pModule;
};

#endif    // MLOPTIMIZER_H
