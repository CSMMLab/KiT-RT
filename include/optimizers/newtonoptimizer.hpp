/*!
 * @file newtonoptimizer.h
 * @brief class for solving the minimal entropy optimization problem using a newton optimizer with line search.
 * @author S. Schotthöfer
 */

#ifndef NEWTONOPTIMIZER_H
#define NEWTONOPTIMIZER_H

#include "optimizerbase.hpp"

class QuadratureBase;

class NewtonOptimizer : public OptimizerBase
{
  public:
    NewtonOptimizer( Config* settings );

    ~NewtonOptimizer();

    void Solve( Vector& alpha, Vector& sol, const VectorVector& moments, unsigned idx_cell = 0 ) override;
    void SolveMultiCell( VectorVector& alpha, VectorVector& sol, const VectorVector& moments ) override;

    /*! @brief Computes the objective function
                grad = <eta(alpha*m)> - alpha*sol */
    virtual double ComputeObjFunc( Vector& alpha, Vector& sol, const VectorVector& moments );

    /*! @brief Computes hessian of objective function and stores it in hessian
        grad = <mXm*eta*'(alpha*m)> */
    virtual void ComputeHessian( Vector& alpha, const VectorVector& moments, Matrix& hessian );

    /*! @brief In 1D, this function scales the quadrature weigths to compute the entropy integrals in arbitrary (bounded) intervals
        @param velocityScale :scaling factor of the symmetric velocity intervall with mean 0
        */
    void ScaleQuadWeights( double velocityScale );

  protected:
    /*! @brief Computes gradient of objective function and stores it in grad
                grad = <m*eta*'(alpha*m)> - sol */
    virtual void ComputeGradient( Vector& alpha, Vector& sol, const VectorVector& moments, Vector& grad );

    QuadratureBase* _quadrature; /*!< @brief used quadrature */    // THis is memory doubling! Try to use a pointer.

    unsigned _nq;    /*!< @brief number of quadrature points */
    Vector _weights; /*!<  @brief quadrature weights, dim(_weights) = (_nq) */

    double _epsilon;                 /*!< @brief Termination criterion for newton optimizer */
    unsigned short _maxIterations;   /*!< @brief Max iterations of the newton solver */
    double _alpha;                   /*!< @brief Newton Step Size */
    unsigned short _maxLineSearches; /*!< @brief Max amount of line searches for Newton Algo */
};

#endif    // NEWTONOPTIMIZER_H
