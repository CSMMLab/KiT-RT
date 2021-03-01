/*!
 * @file optimizerbase.h
 * @brief Base class for solving the minimal entropy optimization problem
 * @author S. Schotthöfer
 */

#ifndef OPTIMIZERBASE_H
#define OPTIMIZERBASE_H

#include "common/globalconstants.h"    // for ML, NEWTON
#include "common/typedef.h"

// Foward declaration
class Config;
class EntropyBase;

class OptimizerBase
{
  public:
    OptimizerBase( Config* settings );

    virtual ~OptimizerBase();

    /*! @brief Optimizer creator: Depending on the chosen option, this function creates an object of the chosen child class of OptimizerBase */
    static OptimizerBase* Create( Config* settings );

    /*! @brief  Computes the optimal Lagrange multilpiers for the dual entropy minimization problem
     *  @param  alpha vector where the solution Lagrange multipliers are saved to.
     *  @param  u  moment vector
     *  @param  moments  VectorVector to the moment basis evaluated at all quadpoints
     *  @param idx_cell index of the cell where alpha should be computed (out of u) */
    virtual void Solve( Vector& alpha, Vector& u, const VectorVector& moments, unsigned idx_cell = 0 ) = 0;

    /*! @brief   Computes the optimal Lagrange multilpiers for the dual entropy minimization problem
     *  @param   alpha vector where the solution Lagrange multipliers are saved to.
     *  @param   u  moment vector
     *  @param   moments  VectorVector to the moment basis evaluated at all quadpoints    */
    virtual void SolveMultiCell( VectorVector& alpha, VectorVector& u, const VectorVector& moments ) = 0;

  protected:
    EntropyBase* _entropy; /*!< @brief Class to handle entropy functional evaluations */
    Config* _settings;     /*!< @brief Pointer to settings class of the solver */
};

#endif    // OPTIMIZERBASE_H
