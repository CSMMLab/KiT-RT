#ifndef OPTIMIZERBASE_H
#define OPTIMIZERBASE_H

#include "entropies/entropybase.h"
#include "settings/config.h"
#include "settings/typedef.h"

class OptimizerBase
{
  public:
    OptimizerBase( Config* settings );

    static OptimizerBase* Create( Config* settings );

    /*! @brief  : Computes the optimal Lagrange multilpiers for the dual entropy minimization problem
     *  @param  : Vector  u = pointer to vector of given moments. // Maybe use pointer for performance?
     *  @return : Vector  alpha = optimal lagrange multipliers. Has the same length as Vector u. */
    virtual Vector Solve( Vector u ) = 0;

  protected:
    EntropyBase* _entropy; /*! @brief: Class to handle entropy functional evaluations */
    Config* _settings;
};

#endif    // OPTIMIZERBASE_H
