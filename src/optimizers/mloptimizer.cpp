
#ifdef BUILD_ML

#include "cppflow/cppflow.h"

#endif

#include "common/config.hpp"
#include "optimizers/mloptimizer.hpp"
#include "toolboxes/errormessages.hpp"
#include <iostream>

MLOptimizer::MLOptimizer( Config* settings ) : OptimizerBase( settings ) {

#ifdef BUILD_ML
    auto input_1 = cppflow::fill( { 10, 2 }, 0.5f );
    auto input_2 = cppflow::fill( { 10, 2 }, 0.5f );
    cppflow::model model( "../ext/neuralEntropy/final_model" );

    auto output =
        model( { { "serving_default_x:0", input_1 } }, { "StatefulPartitionedCall:0", "StatefulPartitionedCall:1", "StatefulPartitionedCall:2" } );

    std::cout << "output_1: " << output[0] << std::endl;
    std::cout << "output_2: " << output[1] << std::endl;
    std::cout << "output_3: " << output[2] << std::endl;
    ErrorMessages::Error( "ML build not configured. Please activate cmake flage BUILD_ML.", CURRENT_FUNCTION );

#else
    ErrorMessages::Error( "ML build not configured. Please activate cmake flage BUILD_ML.", CURRENT_FUNCTION );
#endif
}

MLOptimizer::~MLOptimizer() {}

void MLOptimizer::Solve( Vector& alpha, Vector& u, const VectorVector& /*moments*/, unsigned /*idx_cell*/ ) {}

void MLOptimizer::SolveMultiCell( VectorVector& alpha, VectorVector& u, const VectorVector& /*moments*/ ) {}

void MLOptimizer::ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) {
    ErrorMessages::Error( "This function is not yet implemented.", CURRENT_FUNCTION );
}
