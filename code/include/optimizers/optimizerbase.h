/*!
 * @file optimizerbase.h
 * @brief Base class for solving the minimal entropy optimization problem
 * @author S. Schotthöfer
 */

#ifndef OPTIMIZERBASE_H
#define OPTIMIZERBASE_H

#include "common/typedef.h"

// Foward declaration
class Config;
class EntropyBase;

class OptimizerBase
{
  public:
    OptimizerBase( Config* settings );

    virtual ~OptimizerBase();

    /*! @brief: Optimizer creator: Depending on the chosen option, this function creates an object of the chosen child class of OptimizerBase */
    static OptimizerBase* Create( Config* settings );

    /*! @brief  : Computes the optimal Lagrange multilpiers for the dual entropy minimization problem
     *  @param  : Vector  u = pointer to vector of given moments. // Maybe use pointer for performance?
     *  @return : Vector  alpha = optimal lagrange multipliers. Has the same length as Vector u. */
    virtual void Solve( Vector& lambda, Vector& u, const VectorVector& moments, unsigned idx_cell = 0 ) = 0;

    virtual void SolveMultiCell( VectorVector& lambda, VectorVector& u, const VectorVector& moments ) = 0;

  protected:
    EntropyBase* _entropy; /*! @brief: Class to handle entropy functional evaluations */
    Config* _settings;     /*! @biref: Pointer to settings class of the solver */
};

#endif    // OPTIMIZERBASE_H
