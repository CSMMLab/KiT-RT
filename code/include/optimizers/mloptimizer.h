#ifndef MLOPTIMIZER_H
#define MLOPTIMIZER_H

#include "optimizerbase.h"

class MLOptimizer : public OptimizerBase
{
  public:
    MLOptimizer( Config* settings );

    inline ~MLOptimizer();

    void Solve( Vector& lambda, Vector& u, VectorVector& moments, unsigned idx_cell = 0 ) override;

  private:
    double* callNetwork(const unsigned input_size, double* nn_input);
    void finalize_python();
    void initialize_python();
    void init_numpy();
   // void initialize_Network();

};

#endif    // MLOPTIMIZER_H
