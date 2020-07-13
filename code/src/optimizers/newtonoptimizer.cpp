
#include "optimizers/newtonoptimizer.h"
#include "quadratures/quadraturebase.h"
#include "settings/config.h"

NewtonOptimizer::NewtonOptimizer( Config* settings ) : OptimizerBase( settings ) {
    _quadrature       = QuadratureBase::CreateQuadrature( settings->GetQuadName(), settings->GetQuadOrder() );
    _nq               = _quadrature->GetNq();
    _weights          = _quadrature->GetWeights();
    _quadPointsSphere = _quadrature->GetPointsSphere();
    _maxIterations    = settings->GetMaxIterNewtonOptimizer();
    _alpha            = 0.1;     // settings->GetNewtonStepSize();
    _maxLineSearches  = 1000;    // settings->GetMaxLineSearches();
}

void NewtonOptimizer::ComputeGradient( Vector& alpha, Vector& sol, VectorVector& moments, Vector& grad ) {

    // Reset Vector
    for( unsigned idx_sys; idx_sys < grad.size(); idx_sys++ ) {
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
    for( unsigned idx_Row; idx_Row < moments.size(); idx_Row++ ) {
        for( unsigned idx_Col; idx_Col < moments.size(); idx_Col++ ) {
            hessian( idx_Row, idx_Col ) = 0.0;
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

    // if we have quadratic entropy, then alpha = u;
    if( _settings->GetEntropyName() == QUADRATIC ) {
        lambda = sol;
        return;
    }

    // Start Newton Algorithm
    //
    unsigned nSize = lambda.size();

    // int maxRefinements = 1000;
    // unsigned nTotal    = _nTotalForRef[refLevel];

    Vector grad( nSize, 0.0 );

    // check if initial guess is good enough
    ComputeGradient( lambda, sol, moments, grad );

    // Gradient( g, lambda, u, refLevel );

    if( norm( grad ) < _settings->GetNewtonOptimizerEpsilon() ) {
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

        int lineSearchCounter = 0;

        // Line Search
        while( norm( dlambda ) < norm( dlambdaNew ) || !std::isfinite( norm( dlambdaNew ) ) ) {
            stepSize *= 0.5;

            lambdaNew = lambda - stepSize * _alpha * H * grad;
            ComputeGradient( lambdaNew, sol, moments, dlambdaNew );

            // Check if FONC is fullfilled
            if( norm( dlambdaNew ) < _settings->GetNewtonOptimizerEpsilon() ) {
                lambda = lambdaNew;
                return;
            }
            else if( ++lineSearchCounter > _maxLineSearches ) {
                //_log->error( "[closure] Newton needed too many refinement steps!" );
                std::cout << "[closure] Newton needed too many refinement steps!";
                exit( EXIT_FAILURE );
            }
        }
        lambda = lambdaNew;
        if( norm( dlambdaNew ) < _settings->GetNewtonOptimizerEpsilon() ) {
            lambda = lambdaNew;
            return;
        }
    }
    //    _log->error( "[closure] Newton did not converge!" );
    std::cout << "[closure] Newton did not converge!";
    exit( EXIT_FAILURE );
}
