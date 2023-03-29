#ifndef MNSOLVER_NORMALIZED_H
#define MNSOLVER_NORMALIZED_H

#include "mnsolver.hpp"

class EntropyBase;
class SphericalBase;
class OptimizerBase;

class MNSolverNormalized : public MNSolver
{
  public:
    /**
     * @brief MNSolverNormalized constructor
     * @param settings Config class that stores all needed information
     */
    MNSolverNormalized( Config* settings );

    /*! @brief MNSolverNormalized destructor */
    virtual ~MNSolverNormalized();

  private:
    Vector _u0; /*!< @brief Vector of zeroOrderMoments*/

    // ---- Private Member functions ---

    // Solver
    void IterPreprocessing( unsigned /*idx_iter*/ ) override;

    // Debugging purposes
    OptimizerBase* _optimizer2;
};
#endif    // MNSOLVER_NORMALIZED_H
