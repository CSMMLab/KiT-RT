#ifndef MLOPTIMIZER_H
#define MLOPTIMIZER_H

#include "optimizerbase.hpp"

#ifdef BUILD_ML
#include "cppflow/cppflow.h"

class QuadratureBase;

class NeuralNetworkOptimizer : public OptimizerBase
{
  public:
    NeuralNetworkOptimizer( Config* settings );

    inline ~NeuralNetworkOptimizer();

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

    std::string _tfModelInputName;               /*!< @brief Name of the tf model input */
    cppflow::model* _tfModel;                    /*!< @brief wrapper object for the compiled tensorflow model*/
    cppflow::tensor _modelInput;                 /*!< @brief model input tensor. dims: _nCellsx_nSys*/
    std::vector<float> _modelServingVectorU;     /*!< @brief model input as a 1D vector. dims: _nCells*(_nSys-1) */
    std::vector<float> _modelServingVectorAlpha; /*!< @brief model output as a 1D vector. dims: _nCells*_nSys */

    // std::vector<cppflow::tensor> _modelOutput; /*!< @brief model input tensor. dims: _nModelOutputx_nCellsx_nSys*/
    unsigned _nSystem;                /*!< @brief  size of the moment system including zero order moment*/
    VectorVector _reducedMomentBasis; /*!< @brief reduced basis functions (excluding order zero) */

    vector<Matrix> _rotationMats;  /*!< @brief vector of Rotation matrices for symmetry enforcing */
    vector<Matrix> _rotationMatsT; /*!< @brief vector of transpose Rotation matrices for symmetry enforcing */

    Matrix CreateRotator( const Vector& uFirstMoment ); /*!< @brief Creates a rotation matrix R using the first moment of a momnet vector */
    Vector RotateM1( Vector& vec, Matrix& R );          /*!< @brief Rotates the M1 part of a 2D moment vector using a rotation matrix R */
    /*!< @brief Rotates the tensorized M2 part of a 2D moment vector using a rotation matrix R */
    Matrix RotateM2( Matrix& vec, Matrix& R, Matrix& Rt );
    /*!< @brief Rotates the tensorized M3 part of a 2D moment vector using a rotation matrix R */
    std::vector<VectorVector> RotateM3( std::vector<VectorVector>& vec, Matrix& R );
};
#else
// Dummy class
class NeuralNetworkOptimizer : public OptimizerBase
{
  public:
    NeuralNetworkOptimizer( Config* settings );

    inline ~NeuralNetworkOptimizer();

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
