#ifndef MLOPTIMIZER_H
#define MLOPTIMIZER_H

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "optimizerbase.hpp"

// Forward declaration

class MLOptimizer : public OptimizerBase
{
  public:
    MLOptimizer( Config* settings );

    inline ~MLOptimizer();

    void Solve( Vector& alpha, Vector& u, const VectorVector& moments, unsigned idx_cell = 0 ) override;
    void SolveMultiCell( VectorVector& alpha, VectorVector& u, const VectorVector& moments ) override;

    /*! @brief Reconstruct the moment sol from the Lagrange multiplier alpha
     *  @param sol moment vector
     *  @param alpha Lagrange multipliers
     *  @param moments Moment basis      */
    void ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) override;

  private:

};

#endif    // MLOPTIMIZER_H
