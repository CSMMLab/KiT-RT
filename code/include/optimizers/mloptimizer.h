#ifndef MLOPTIMIZER_H
#define MLOPTIMIZER_H

#include "optimizerbase.h"

class MLOptimizer : public OptimizerBase
{
  public:
    MLOptimizer( Config* settings );

    inline ~MLOptimizer() {}

    void Solve( Vector& lambda, Vector& u, VectorVector& moments, unsigned idx_cell = 0 ) override;

  private:
};

#endif    // MLOPTIMIZER_H
