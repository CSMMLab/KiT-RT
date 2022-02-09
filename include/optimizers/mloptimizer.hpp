#ifndef MLOPTIMIZER_H
#define MLOPTIMIZER_H

#include "optimizerbase.hpp"

#ifdef BUILD_ML
#include "cppflow/cppflow.h"

class QuadratureBase;

class MLOptimizer : public OptimizerBase
{
  public:
    MLOptimizer( Config* settings );

    inline ~MLOptimizer();

    void Solve( Vector& alpha, Vector& u, const VectorVector& moments, unsigned idx_cell = 0 ) override;

    void SolveMultiCell( VectorVector& alpha, VectorVector& u, const VectorVector& moments ) override;

    /*! @brief Reconstruct the moment sol from the Lagrange multiplier alpha
     *  @param sol moment vector
     *  @param alpha Lagrange multipliers
     *  @param moments Moment basis      */
    void ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) override;

  private:
    QuadratureBase* _quadrature; /*!< @brief used quadrature */    // THis is memory doubling! Try to use a pointer.
    unsigned _nq;                                                  /*!< @brief number of quadrature points */
    Vector _weights;                                               /*!<  @brief quadrature weights, dim(_weights) = (_nq) */

    cppflow::model* _tfModel;                    /*!< @brief wrapper object for the compiled tensorflow model*/
    cppflow::tensor _modelInput;                 /*!< @brief model input tensor. dims: _nCellsx_nSys*/
    std::vector<float> _modelServingVectorU;     /*!< @brief model input as a 1D vector. dims: _nCells*(_nSys-1) */
    std::vector<float> _modelServingVectorAlpha; /*!< @brief model output as a 1D vector. dims: _nCells*_nSys */

    // std::vector<cppflow::tensor> _modelOutput; /*!< @brief model input tensor. dims: _nModelOutputx_nCellsx_nSys*/
    unsigned _nSystem;                /*!< @brief  size of the moment system including zero order moment*/
    VectorVector _reducedMomentBasis; /*!< @brief reduced basis functions (excluding order zero) */
};
#else
// Dummy class
class MLOptimizer : public OptimizerBase
{
  public:
    MLOptimizer( Config* settings );

    inline ~MLOptimizer();

    inline void Solve( Vector& alpha, Vector& u, const VectorVector& moments, unsigned idx_cell = 0 ) override{};

    inline void SolveMultiCell( VectorVector& alpha, VectorVector& u, const VectorVector& moments ) override{};

    /*! @brief Reconstruct the moment sol from the Lagrange multiplier alpha
     *  @param sol moment vector
     *  @param alpha Lagrange multipliers
     *  @param moments Moment basis      */
    inline void ReconstructMoments( Vector& sol, const Vector& alpha, const VectorVector& moments ) override{};
};
#endif
#endif    // MLOPTIMIZER_H
