/*!
 * @file partregularizednewtonoptimizer.h
 * @brief class for solving the minimal entropy optimization problem using a partly regularized (order zero moment is not regularized) newton
 * optimizer with line search.
 * @author S. Schotthöfer
 */

#ifndef PARTREGULARIZEDNEWTONOPTIMIZER_H
#define PARTREGULARIZEDNEWTONOPTIMIZER_H

#include "newtonoptimizer.hpp"

class QuadratureBase;

class PartRegularizedNewtonOptimizer : public NewtonOptimizer
{
  public:
    PartRegularizedNewtonOptimizer( Config* settings );

    ~PartRegularizedNewtonOptimizer();

    /*! @brief Computes the objective function
                grad = <eta(alpha*m)> - alpha*sol  + gamma/2*norm(alpha)*/
    double ComputeObjFunc( Vector& alpha, Vector& sol, const VectorVector& moments ) override;

    /*! @brief Computes hessian of objective function and stores it in hessian
        grad = <mXm*eta*'(alpha*m)> */
    void ComputeHessian( Vector& alpha, const VectorVector& moments, Matrix& hessian ) override;

    /*! @brief In 1D, this function scales the quadrature weigths to compute the entropy integrals in arbitrary (bounded) intervals
        @param leftBound : left boundary of the interval
        @param rightBound : right boundary of the interval*/
    void ScaleQuadWeights( double leftBound, double rightBound );

  private:
    /*! @brief Computes gradient of objective function and stores it in grad
                grad = <m*eta*'(alpha*m)> - sol + _gamma*alpha */
    void ComputeGradient( Vector& alpha, Vector& sol, const VectorVector& moments, Vector& grad ) override;

    double _gamma; /*!<  @brief Regularization Parameter*/
};

#endif    // PARTREGULARIZEDNEWTONOPTIMIZER_H
