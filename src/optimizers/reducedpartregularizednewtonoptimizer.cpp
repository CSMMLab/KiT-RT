/*!
 * @file reducedpartregularizednewtonoptimizer.cpp
 * @brief class for solving the minimal entropy optimization problem using a partly regularized (order zero moment is not regularized) newton
 * optimizer with line search.
 * @author S. Schotthöfer
 */

#include "optimizers/reducedpartregularizednewtonoptimizer.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

ReducedPartRegularizedNewtonOptimizer::ReducedPartRegularizedNewtonOptimizer( Config* settings ) : ReducedNewtonOptimizer( settings ) {
    _gamma = settings->GetRegularizerGamma();    // Regularization parameter (to be included in settings)
}

ReducedPartRegularizedNewtonOptimizer::~ReducedPartRegularizedNewtonOptimizer() {}

double ReducedPartRegularizedNewtonOptimizer::ComputeObjFunc( Vector& alpha, Vector& sol, const VectorVector& moments ) {
    double result = ReducedNewtonOptimizer::ComputeObjFunc( alpha, sol, moments );    // Calls non regularized objective function
    for( unsigned idx_sys = 0; idx_sys < alpha.size(); idx_sys++ ) {
        result += 0.5 * _gamma * alpha[idx_sys] * alpha[idx_sys];    // Add regularizer norm(_alpha^r)^2
    }
    return result;
}

void ReducedPartRegularizedNewtonOptimizer::ComputeGradient( Vector& alpha, Vector& sol, const VectorVector& moments, Vector& grad ) {
    ReducedNewtonOptimizer::ComputeGradient( alpha, sol, moments, grad );    // compute unregularized gradients
    for( unsigned idx_sys = 0; idx_sys < alpha.size(); idx_sys++ ) {
        grad[idx_sys] += _gamma * alpha[idx_sys];    // Add regularizer _alpha^r
    }
}

void ReducedPartRegularizedNewtonOptimizer::ComputeHessian( Vector& alpha, const VectorVector& moments, Matrix& hessian ) {
    ReducedNewtonOptimizer::ComputeHessian( alpha, moments, hessian );    // compute unregularized hessian)
    for( unsigned idx_sys = 0; idx_sys < alpha.size(); idx_sys++ ) {
        hessian( idx_sys, idx_sys ) += _gamma;    // Add block identity matrix with regularizer
    }
}

void ReducedPartRegularizedNewtonOptimizer::ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) {
    ReducedNewtonOptimizer::ReconstructMoments( sol, alpha, moments );
    for( unsigned idx_sys = 0; idx_sys < alpha.size(); idx_sys++ ) {
        sol[idx_sys] += _gamma * alpha[idx_sys];    // Add regularizer _alpha^r
    }
}
