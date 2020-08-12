#ifndef OPTIMIZERBASE_H
#define OPTIMIZERBASE_H

#include "common/typedef.h"
#include "entropies/entropybase.h"

// Foward declaration
class Config;

class OptimizerBase
{
  public:
    OptimizerBase( Config* settings );

    virtual inline ~OptimizerBase() { delete _entropy; }

    static OptimizerBase* Create( Config* settings );

    /*! @brief  : Computes the optimal Lagrange multilpiers for the dual entropy minimization problem
     *  @param  : Vector  u = pointer to vector of given moments. // Maybe use pointer for performance?
     *  @return : Vector  alpha = optimal lagrange multipliers. Has the same length as Vector u. */
    virtual void Solve( Vector& lambda, Vector& u, VectorVector& moments, unsigned idx_cell = 0 ) = 0;

    virtual void SolveMultiCell( VectorVector& lambda, VectorVector& u, VectorVector& moments ) = 0;

  protected:
    EntropyBase* _entropy; /*! @brief: Class to handle entropy functional evaluations */
    Config* _settings;
};

#endif    // OPTIMIZERBASE_H
