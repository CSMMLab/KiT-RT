/*!
 * @file newtonoptimizer.cpp
 * @brief class for solving the minimal entropy optimization problem using a newton optimizer with line search.
 * @author S. Schotthöfer
 */

#include "optimizers/newtonoptimizer.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

#include <omp.h>

NewtonOptimizer::NewtonOptimizer( Config* settings ) : OptimizerBase( settings ) {
    _quadrature      = QuadratureBase::Create( settings );
    _nq              = _quadrature->GetNq();
    _weights         = _quadrature->GetWeights();
    _maxIterations   = settings->GetNewtonIter();
    _delta           = settings->GetNewtonStepSize();
    _maxLineSearches = settings->GetNewtonMaxLineSearches();
    _epsilon         = settings->GetNewtonOptimizerEpsilon();
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

void NewtonOptimizer::SolveMultiCell( VectorVector& alpha, VectorVector& sol, const VectorVector& moments, Vector& alpha_norms ) {

    unsigned nCells = alpha.size();
    unsigned nSys = alpha[0].size();

    // if we  have quadratic entropy, then alpha = u;
    if( _settings->GetEntropyName() == QUADRATIC && _settings->GetNewtonFastMode() ) {
        for( unsigned idx_cell = 0; idx_cell < nCells; idx_cell++ ) {
            alpha[idx_cell] = sol[idx_cell];
        }
        return;
    }

#pragma omp parallel for schedule( guided )
    for( unsigned idx_cell = 0; idx_cell < nCells; idx_cell++ ) {
        Solve( alpha[idx_cell], sol[idx_cell], moments, idx_cell );
        alpha_norms[idx_cell] = 0.0;
        for ( unsigned idx_sys = 1; idx_sys < nSys; idx_sys++ ){
            alpha_norms[idx_cell] += alpha[idx_cell][idx_sys]*alpha[idx_cell][idx_sys];
        }
    }
}

void NewtonOptimizer::Solve( Vector& alpha, Vector& sol, const VectorVector& moments, unsigned idx_cell ) {

    /* solve the problem argmin ( <eta_*(alpha*m)>-alpha*u))
     * where alpha = Lagrange multiplier
     *           m = moment basis
     *           u = current "moment solution"
     */

    // if we  have quadratic entropy, then alpha = u;
    if( _settings->GetEntropyName() == QUADRATIC && _settings->GetNewtonFastMode() ) {
        alpha = sol;
        return;
    }

    // Start Newton Algorithm

    unsigned nSize = alpha.size();

    Vector grad( nSize, 0.0 );

    // check if initial guess is good enough
    ComputeGradient( alpha, sol, moments, grad );

    if( norm( grad ) < _epsilon ) {
        return;
    }

    // If not, compute the Hessian
    Matrix H( nSize, nSize, 0.0 );

    // calculate initial Hessian and gradient
    Vector dalpha = -grad;
    ComputeHessian( alpha, moments, H );

    // invert Hessian
    invert( H );

    if( _maxIterations == 1 ) {
        alpha = alpha - _delta * H * grad;
        return;
    }

    Vector alphaNew( nSize, 0.0 );    //  New newton step
    Vector dalphaNew( nSize );        //  Gradient at New newton step

    alphaNew = alpha - _delta * H * grad;

    // Compute Gradient of new point;
    ComputeGradient( alphaNew, sol, moments, dalphaNew );

    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            ComputeGradient( alpha, sol, moments, grad );

            dalpha = -grad;
            ComputeHessian( alpha, moments, H );
            invert( H );
            alphaNew = alpha - _delta * H * grad;
            ComputeGradient( alphaNew, sol, moments, dalphaNew );
        }

        // Line Search

        int lineSearchCounter = 0;

        while( norm( dalpha ) < norm( dalphaNew ) || !std::isfinite( norm( dalphaNew ) ) ) {
            stepSize *= 0.5;

            alphaNew = alpha - stepSize * _delta * H * grad;
            ComputeGradient( alphaNew, sol, moments, dalphaNew );

            // Check if FONC is locally fullfilled
            if( norm( dalphaNew ) < _epsilon ) {
                alpha = alphaNew;
                return;
            }
            else if( ++lineSearchCounter > _maxLineSearches ) {
                std::string uSolString = "At moment: (" + std::to_string( sol[0] );
                for( unsigned i = 1; i < nSize; i++ ) {
                    uSolString += " | " + std::to_string( sol[i] );
                }
                uSolString += ").";
                std::string alphaSolString = "At mulitplier: (" + std::to_string( alpha[0] );
                for( unsigned i = 1; i < nSize; i++ ) {
                    alphaSolString += " | " + std::to_string( alpha[i] );
                }
                alphaSolString += ").";
                ErrorMessages::Error( "Newton needed too many refinement steps!  at cell " + std::to_string( idx_cell ) + ",\n " + uSolString + "\n" +
                                          alphaSolString,
                                      CURRENT_FUNCTION );
            }
        }
        alpha = alphaNew;
        if( norm( dalphaNew ) < _epsilon ) {
            alpha = alphaNew;
            return;
        }
    }
    std::string uSolString = "At moment: (" + std::to_string( sol[0] );
    for( unsigned i = 1; i < nSize; i++ ) {
        uSolString += " | " + std::to_string( sol[i] );
    }
    uSolString += ").";

    std::string alphaSolString = "At mulitplier: (" + std::to_string( alpha[0] );
    for( unsigned i = 1; i < nSize; i++ ) {
        alphaSolString += " | " + std::to_string( alpha[i] );
    }
    alphaSolString += ").";

    if( _settings->GetDim() == 3 ) {
        Vector u1     = { sol[1], sol[2], sol[3] };
        double normU1 = norm( u1 );
        ErrorMessages::Error( "Newton did not converge at cell " + std::to_string( idx_cell ) + "\n" + uSolString + "\n " + alphaSolString +
                                  "\nObjective function value: " + std::to_string( ComputeObjFunc( alpha, sol, moments ) ) +
                                  "\nBoundary Ratio: " + std::to_string( normU1 / sol[0] ),
                              CURRENT_FUNCTION );
    }
    else {
        ErrorMessages::Error( "Newton did not converge at cell " + std::to_string( idx_cell ) + "\n" + uSolString + "\nMultiplier: " +
                                  alphaSolString + "\nObjective function value: " + std::to_string( ComputeObjFunc( alpha, sol, moments ) ),
                              CURRENT_FUNCTION );
    }
}

void NewtonOptimizer::ScaleQuadWeights( double velocityScale ) {
    _quadrature->ScalePointsAndWeights( velocityScale );
    _weights = _quadrature->GetWeights();
}

void NewtonOptimizer::ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) {
    double entropyReconstruction = 0.0;
    for( unsigned idx_sys = 0; idx_sys < sol.size(); idx_sys++ ) {
        sol[idx_sys] = 0.0;
    }
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        // Make entropyReconstruction a member vector, s.t. it does not have to be re-evaluated in ConstructFlux
        entropyReconstruction = _entropy->EntropyPrimeDual( blaze::dot( alpha, moments[idx_quad] ) );
        sol += moments[idx_quad] * ( _weights[idx_quad] * entropyReconstruction );
    }
}
