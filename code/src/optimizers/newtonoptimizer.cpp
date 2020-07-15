
#include "optimizers/newtonoptimizer.h"
#include "common/config.h"
#include "quadratures/quadraturebase.h"
#include "toolboxes/errormessages.h"

NewtonOptimizer::NewtonOptimizer( Config* settings ) : OptimizerBase( settings ) {
    _quadrature       = QuadratureBase::CreateQuadrature( settings->GetQuadName(), settings->GetQuadOrder() );
    _nq               = _quadrature->GetNq();
    _weights          = _quadrature->GetWeights();
    _quadPointsSphere = _quadrature->GetPointsSphere();
    _maxIterations    = settings->GetNewtonIter();
    _alpha            = settings->GetNewtonStepSize();
    _maxLineSearches  = settings->GetMaxLineSearches();
    _epsilon          = settings->GetNewtonOptimizerEpsilon();
}

double NewtonOptimizer::ComputeObjFunc( Vector& alpha, Vector& sol, VectorVector& moments ) {
    double result = 0.0;

    // Integrate
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        result += _entropy->EntropyDual( dot( alpha, moments[idx_quad] ) ) * _weights[idx_quad];
    }
    result -= dot( alpha, sol );
    return result;
}

void NewtonOptimizer::ComputeGradient( Vector& alpha, Vector& sol, VectorVector& moments, Vector& grad ) {

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

void NewtonOptimizer::ComputeHessian( Vector& alpha, VectorVector& moments, Matrix& hessian ) {
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

void NewtonOptimizer::Solve( Vector& lambda, Vector& sol, VectorVector& moments ) {

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
                ErrorMessages::Error( "Newton needed too many refinement steps!", CURRENT_FUNCTION );
            }
        }
        lambda = lambdaNew;
        if( norm( dlambdaNew ) < _epsilon ) {
            lambda = lambdaNew;
            return;
        }
    }
    ErrorMessages::Error( "Newton did not converge!", CURRENT_FUNCTION );
}
