/*!
 * @file newtonoptimizer.h
 * @brief class for solving the minimal entropy optimization problem using a newton optimizer with line search.
 * @author S. Schotthöfer
 */

#ifndef REDUCEDNEWTONOPTIMIZER_H
#define REDUCEDNEWTONOPTIMIZER_H

#include "newtonoptimizer.hpp"

class ReducedNewtonOptimizer : public NewtonOptimizer
{
  public:
    ReducedNewtonOptimizer( Config* settings );

    ~ReducedNewtonOptimizer();

    /*! @brief Computes the objective function
                grad = <eta(alpha*m)> - alpha*sol */
    virtual double ComputeObjFunc( const Vector& alpha, const Vector& sol, const VectorVector& moments ) override;

    /*! @brief Computes hessian of objective function and stores it in hessian
        grad = <mXm*eta*'(alpha*m)> */
    virtual void ComputeHessian( const Vector& alpha, const VectorVector& moments, Matrix& hessian ) override;

    /*! @brief Reconstruct the moment sol from the Lagrange multiplier alpha
     *  @param sol moment vector
     *  @param alpha Lagrange multipliers
     *  @param moments Moment basis
     */
    virtual void ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) override;

  protected:
    /*! @brief Computes gradient of objective function and stores it in grad
                grad = <m*eta*'(alpha*m)> - sol */
    virtual void ComputeGradient( const Vector& alpha, const Vector& sol, const VectorVector& moments, Vector& grad ) override;
};

#endif    // REDUCEDNEWTONOPTIMIZER_H
