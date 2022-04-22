/*!
 * @file newtonoptimizer.cpp
 * @brief class for solving the minimal entropy optimization problem using a regularized newton optimizer with line search.
 * @author S. Schotthöfer
 */

#include "optimizers/regularizednewtonoptimizer.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "optimizers/newtonoptimizer.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

#include <omp.h>

RegularizedNewtonOptimizer::RegularizedNewtonOptimizer( Config* settings ) : NewtonOptimizer( settings ) {
    _gamma = settings->GetRegularizerGamma();    // Regularization parameter (to be included in settings)
}

RegularizedNewtonOptimizer::~RegularizedNewtonOptimizer() {}

double RegularizedNewtonOptimizer::ComputeObjFunc( Vector& alpha, Vector& sol, const VectorVector& moments ) {
    double result = NewtonOptimizer::ComputeObjFunc( alpha, sol, moments );    // Calls non regularized objective function
    result += 0.5 * _gamma * dot( alpha, alpha );                              // Add regularizer norm(_alpha)
    return result;
}

void RegularizedNewtonOptimizer::ComputeGradient( Vector& alpha, Vector& sol, const VectorVector& moments, Vector& grad ) {
    NewtonOptimizer::ComputeGradient( alpha, sol, moments, grad );    // compute unregularized gradients
    grad += _gamma * alpha;
}

void RegularizedNewtonOptimizer::ComputeHessian( Vector& alpha, const VectorVector& moments, Matrix& hessian ) {
    NewtonOptimizer::ComputeHessian( alpha, moments, hessian );    // compute unregularized hessian)
    hessian += _gamma * IdentityMatrix( alpha.size() );
}

void RegularizedNewtonOptimizer::ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) {
    NewtonOptimizer::ReconstructMoments( sol, alpha, moments );
    sol += _gamma * alpha;    // Add regularizer _alpha^r
}
