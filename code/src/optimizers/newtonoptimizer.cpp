
#include "optimizers/newtonoptimizer.h"
#include "settings/config.h"

void NewtonOptimizer::Solve( Vector& lambda, Vector& u ) {

    // if we have quadratic entropy, then alpha = u;
    if( _settings->GetEntropyName() == QUADRATIC ) return;

    // Start Newton Algorithm
    //
    //    int maxRefinements = 1000;
    //    unsigned nTotal    = _nTotalForRef[refLevel];
    //
    //    Vector g( _nStates * nTotal );
    //
    //    // check if initial guess is good enough
    //    Gradient( g, lambda, u, refLevel );
    //    if( CalcNorm( g, nTotal ) < _settings->GetEpsilon() ) {
    //        return;
    //    }
    //    Matrix H( _nStates * nTotal, _nStates * nTotal );
    //    Vector dlambdaNew( _nStates * nTotal );
    //    // calculate initial Hessian and gradient
    //    Vector dlambda = -g;
    //    // std::cout << g << std::endl;
    //    Hessian( H, lambda, refLevel );
    //
    //    // double eps     = 0.00001;
    //    // lambda( 1, 0 ) = lambda( 1, 0 ) + eps;
    //
    //    // Vector gEps( _nStates * nTotal );
    //
    //    // Gradient( gEps, lambda, u, refLevel );
    //
    //    // std::cout << "H_FD = " << ( gEps - g ) / eps << std::endl;
    //    // std::cout << "H = " << H << std::endl;
    //    // exit( EXIT_FAILURE );
    //
    //    posv( H, g );
    //    if( _maxIterations == 1 ) {
    //        AddMatrixVectorToMatrix( lambda, -_alpha * g, lambda, nTotal );
    //        return;
    //    }
    //    Matrix lambdaNew( _nStates, nTotal );
    //    AddMatrixVectorToMatrix( lambda, -_alpha * g, lambdaNew, nTotal );
    //    Gradient( dlambdaNew, lambdaNew, u, refLevel );
    //    // perform Newton iterations
    //    for( unsigned l = 0; l < _maxIterations; ++l ) {
    //        double stepSize = 1.0;
    //        if( l != 0 ) {
    //            Gradient( g, lambda, u, refLevel );
    //            dlambda = -g;
    //            Hessian( H, lambda, refLevel );
    //            posv( H, g );
    //            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
    //            Gradient( dlambdaNew, lambdaNew, u, refLevel );
    //        }
    //        int refinementCounter = 0;
    //        // std::cout << CalcNorm( dlambdaNew, nTotal ) << std::endl;
    //        while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) || !std::isfinite( CalcNorm( dlambdaNew, nTotal ) ) ) {
    //            stepSize *= 0.5;
    //            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
    //            Gradient( dlambdaNew, lambdaNew, u, refLevel );
    //            if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
    //                lambda = lambdaNew;
    //                return;
    //            }
    //            else if( ++refinementCounter > maxRefinements ) {
    //                _log->error( "[closure] Newton needed too many refinement steps!" );
    //                exit( EXIT_FAILURE );
    //            }
    //        }
    //        lambda = lambdaNew;
    //        if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
    //            lambda = lambdaNew;
    //            return;
    //        }
    //    }
    //    _log->error( "[closure] Newton did not converge!" );
    //    exit( EXIT_FAILURE );

    // return Vector( 1, 0.0 );    // dummy
}
