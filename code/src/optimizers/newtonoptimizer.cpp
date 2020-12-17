/*!
 * @file newtonoptimizer.cpp
 * @brief class for solving the minimal entropy optimization problem using a newton optimizer with line search.
 * @author S. Schotthöfer
 */

#include "optimizers/newtonoptimizer.h"
#include "common/config.h"
#include "entropies/entropybase.h"
#include "quadratures/quadraturebase.h"
#include "toolboxes/errormessages.h"

#include <omp.h>

NewtonOptimizer::NewtonOptimizer( Config* settings ) : OptimizerBase( settings ) {
    _quadrature       = QuadratureBase::Create( settings );
    _nq               = _quadrature->GetNq();
    _weights          = _quadrature->GetWeights();
    _quadPointsSphere = _quadrature->GetPointsSphere();
    _maxIterations    = settings->GetNewtonIter();
    _alpha            = settings->GetNewtonStepSize();
    _maxLineSearches  = settings->GetNewtonMaxLineSearches();
    _epsilon          = settings->GetNewtonOptimizerEpsilon();
}

NewtonOptimizer::~NewtonOptimizer() { delete _quadrature; }

double NewtonOptimizer::ComputeObjFunc( Vector& alpha, Vector& sol, const VectorVector& moments ) {
    double result = 0.0;

    // Integrate
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        result += _entropy->EntropyDual( dot( alpha, moments[idx_quad] ) ) * _weights[idx_quad];
    }
    result -= dot( alpha, sol );
    return result;
}

void NewtonOptimizer::ComputeGradient( Vector& alpha, Vector& sol, const VectorVector& moments, Vector& grad ) {

    // Reset Vector
    for( unsigned idx_sys = 0; idx_sys < grad.size(); idx_sys++ ) {
        grad[idx_sys] = 0.0;
    }

    // Integrate
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        grad += moments[idx_quad] * ( _entropy->EntropyPrimeDual( dot( alpha, moments[idx_quad] ) ) * _weights[idx_quad] );
    }
    grad -= sol;
}

void NewtonOptimizer::ComputeHessian( Vector& alpha, const VectorVector& moments, Matrix& hessian ) {
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
}

void NewtonOptimizer::SolveMultiCell( VectorVector& lambda, VectorVector& sol, const VectorVector& moments ) {

    unsigned nCells = lambda.size();

    // if we  have quadratic entropy, then alpha = u;
    if( _settings->GetEntropyName() == QUADRATIC && _settings->GetNewtonFastMode() ) {
        for( unsigned idx_cell = 0; idx_cell < nCells; idx_cell++ ) {
            lambda[idx_cell] = sol[idx_cell];
        }
        return;
    }

#pragma omp parallel for schedule( guided )
    for( unsigned idx_cell = 0; idx_cell < nCells; idx_cell++ ) {
        Solve( lambda[idx_cell], sol[idx_cell], moments, idx_cell );
    }
}

void NewtonOptimizer::Solve( Vector& lambda, Vector& sol, const VectorVector& moments, unsigned idx_cell ) {

    /* solve the problem argmin ( <eta(alpha*m)>-alpha*u))
     * where alpha = Lagrange multiplier
     *           m = moment basis
     *           u = current "moment solution"
     */

    // if we  have quadratic entropy, then alpha = u;
    if( _settings->GetEntropyName() == QUADRATIC && _settings->GetNewtonFastMode() ) {
        lambda = sol;
        return;
    }

    // Start Newton Algorithm

    unsigned nSize = lambda.size();

    Vector grad( nSize, 0.0 );

    // check if initial guess is good enough
    ComputeGradient( lambda, sol, moments, grad );

    if( norm( grad ) < _epsilon ) {
        return;
    }

    // If not, compute the Hessian
    Matrix H( nSize, nSize, 0.0 );

    // calculate initial Hessian and gradient
    Vector dlambda = -grad;
    ComputeHessian( lambda, moments, H );

    // invert Hessian
    invert( H );

    if( _maxIterations == 1 ) {
        lambda = lambda - _alpha * H * grad;
        return;
    }

    Vector lambdaNew( nSize, 0.0 );    //  New newton step
    Vector dlambdaNew( nSize );        //  Gradient at New newton step

    lambdaNew = lambda - _alpha * H * grad;

    // Compute Gradient of new point;
    ComputeGradient( lambdaNew, sol, moments, dlambdaNew );

    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            ComputeGradient( lambda, sol, moments, grad );

            dlambda = -grad;
            ComputeHessian( lambda, moments, H );
            invert( H );
            lambdaNew = lambda - _alpha * H * grad;
            ComputeGradient( lambdaNew, sol, moments, dlambdaNew );
        }

        // Line Search

        int lineSearchCounter = 0;

        while( norm( dlambda ) < norm( dlambdaNew ) || !std::isfinite( norm( dlambdaNew ) ) ) {
            stepSize *= 0.5;

            lambdaNew = lambda - stepSize * _alpha * H * grad;
            ComputeGradient( lambdaNew, sol, moments, dlambdaNew );

            // Check if FONC is fullfilled
            if( norm( dlambdaNew ) < _epsilon ) {
                lambda = lambdaNew;
                return;
            }
            else if( ++lineSearchCounter > _maxLineSearches ) {
                ErrorMessages::Error( "Newton needed too many refinement steps!  at cell " + std::to_string( idx_cell ), CURRENT_FUNCTION );
            }
        }
        lambda = lambdaNew;
        if( norm( dlambdaNew ) < _epsilon ) {
            lambda = lambdaNew;
            return;
        }
    }
    ErrorMessages::Error( " Newton did not converge! Norm of gradient is: " + std::to_string( norm( dlambdaNew ) ) + " at cell " +
                              std::to_string( idx_cell ) + ".\nObjective function value is " +
                              std::to_string( ComputeObjFunc( lambda, sol, moments ) ) + " .",
                          CURRENT_FUNCTION );
}
