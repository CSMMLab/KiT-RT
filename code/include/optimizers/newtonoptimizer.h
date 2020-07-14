#ifndef NEWTONOPTIMIZER_H
#define NEWTONOPTIMIZER_H

#include "optimizerbase.h"

class QuadratureBase;

class NewtonOptimizer : public OptimizerBase
{
  public:
    NewtonOptimizer( Config* settings );

    inline ~NewtonOptimizer() {}

    void Solve( Vector& lambda, Vector& u, VectorVector& moments ) override;

  private:
    /*! @brief: Computes gradient of objective function and stores it in grad
                grad = <m*eta*'(alpha*m)> - sol */
    void ComputeGradient( Vector& alpha, Vector& sol, VectorVector& moments, Vector& grad );

    /*! @brief: Computes hessian of objective function and stores it in hessian
                grad = <mXm*eta*'(alpha*m)> */
    void ComputeHessian( Vector& alpha, VectorVector& moments, Matrix& hessian );

    double ComputeObjFunc( Vector& alpha, Vector& sol, VectorVector& moments );

    QuadratureBase* _quadrature; /*! @brief: used quadrature */    // THis is memory doubling! Try to use a pointer.
    unsigned _nq;                                                  /*! @brief: number of quadrature points */
    Vector _weights;                                               /*!  @brief quadrature weights, dim(_weights) = (_nq) */
    VectorVector _quadPointsSphere;                                /*!  @brief (my,phi), dim(_quadPoints) = (_nq,2) */

    double _epsilon;                 /*! @brief: Termination criterion for newton optimizer */
    unsigned short _maxIterations;   /*! @brief: Max iterations of the newton solver */
    double _alpha;                   /*! @brief: Newton Step Size */
    unsigned short _maxLineSearches; /*! @brief: Max amount of line searches for Newton Algo */
};

#endif    // NEWTONOPTIMIZER_H
