/*!
 * @file reducednewtonoptimizer.cpp
 * @brief class for solving the minimal entropy optimization problem using a newton optimizer with line search.
 * @author S. Schotthöfer
 */

#include "optimizers/reducednewtonoptimizer.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"

ReducedNewtonOptimizer::ReducedNewtonOptimizer( Config* settings ) : NewtonOptimizer( settings ) {}

ReducedNewtonOptimizer::~ReducedNewtonOptimizer() { delete _quadrature; }

double ReducedNewtonOptimizer::ComputeObjFunc( const Vector& alpha, const Vector& sol, const VectorVector& moments ) {
    double result = 0.0;

    // Integrate
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        result += _entropy->EntropyDual( dot( alpha, moments[idx_quad] ) ) * _weights[idx_quad];
    }
    // log
    result = log( result );
    // linear term
    result -= dot( alpha, sol );
    return result + 1.0;
}

void ReducedNewtonOptimizer::ComputeGradient( const Vector& alpha, const Vector& sol, const VectorVector& moments, Vector& grad ) {

    // Reset Vector
    for( unsigned idx_sys = 0; idx_sys < grad.size(); idx_sys++ ) {
        grad[idx_sys] = 0.0;
    }

    // Integrate
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        grad += moments[idx_quad] * ( _entropy->EntropyPrimeDual( dot( alpha, moments[idx_quad] ) ) * _weights[idx_quad] );
    }

    double result = 0.0;

    // Integrate
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        result += _entropy->EntropyDual( dot( alpha, moments[idx_quad] ) ) * _weights[idx_quad];
    }

    grad /= result;    // normalization

    // linear part
    grad -= sol;
}

void ReducedNewtonOptimizer::ComputeHessian( const Vector& alpha, const VectorVector& moments, Matrix& hessian ) {
    // Reset Matrix
    unsigned nSize = alpha.size();

    for( unsigned idx_Row = 0; idx_Row < nSize; idx_Row++ ) {
        for( unsigned idx_Col = 0; idx_Col < nSize; idx_Col++ ) {
            hessian( idx_Row, idx_Col ) = 0.0;
            // if( idx_Col == idx_Row ) hessian( idx_Row, idx_Col ) = 1.0;
        }
    }
    // Integrate
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        hessian +=
            outer( moments[idx_quad], moments[idx_quad] ) * ( _entropy->EntropyHessianDual( dot( alpha, moments[idx_quad] ) ) * _weights[idx_quad] );
    }

    double result = 0.0;

    // Integrate
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        result += _entropy->EntropyDual( dot( alpha, moments[idx_quad] ) ) * _weights[idx_quad];
    }

    hessian /= result;    // normalization

    // second part
    Vector grad( nSize, 0.0 );
    // Integrate
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        grad += moments[idx_quad] * ( _entropy->EntropyPrimeDual( dot( alpha, moments[idx_quad] ) ) * _weights[idx_quad] );
    }

    hessian += -1 * outer( grad, grad ) / ( result * result );
}

void ReducedNewtonOptimizer::ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) {
    double entropyReconstruction = 0.0;
    for( unsigned idx_sys = 0; idx_sys < sol.size(); idx_sys++ ) {
        sol[idx_sys] = 0.0;
    }
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        // Make entropyReconstruction a member vector, s.t. it does not have to be re-evaluated in ConstructFlux
        entropyReconstruction = _entropy->EntropyPrimeDual( blaze::dot( alpha, moments[idx_quad] ) );
        sol += moments[idx_quad] * ( _weights[idx_quad] * entropyReconstruction );
    }

    double result = 0.0;

    // Integrate
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        result += _entropy->EntropyDual( dot( alpha, moments[idx_quad] ) ) * _weights[idx_quad];
    }

    sol /= result;    // normalization
}
