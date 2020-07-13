#ifndef NEWTONOPTIMIZER_H
#define NEWTONOPTIMIZER_H

#include "optimizerbase.h"

class NewtonOptimizer : public OptimizerBase
{
  public:
    inline NewtonOptimizer( Config* settings ) : OptimizerBase( settings ) {}

    inline ~NewtonOptimizer() {}

    void Solve( Vector& lambda, Vector& u ) override;
};

#endif    // NEWTONOPTIMIZER_H
