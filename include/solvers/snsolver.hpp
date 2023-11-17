#ifndef SNSOLVER_H
#define SNSOLVER_H

#include "solvers/solverbase.hpp"

class SNSolver : public SolverBase
{
protected:
  Matrix _scatteringKernel; /*!<  @brief scattering kernel for the quadrature */

  // quadrature related numbers

  VectorVector _quadPoints; /*!<  @brief quadrature points, dim(_quadPoints) =
                               (_nq,spatialDim) */
  Vector _weights;          /*!<  @brief quadrature weights, dim(_weights) = (_nq) */

public:
  /*! @brief SNSolver constructor
   * @param settings stores all needed information
   */
  SNSolver(Config *settings);

  virtual ~SNSolver() {}

protected:
  // IO
  void virtual PrepareVolumeOutput() override;
  void virtual WriteVolumeOutput(unsigned idx_pseudoTime) override;

  // Solver
  void virtual FVMUpdate(unsigned idx_energy) override;
  void virtual FluxUpdate() override;
  void virtual IterPreprocessing(unsigned idx_pseudotime) override;
  void virtual IterPostprocessing(unsigned idx_pseudotime) override;

  // Helper
  void ComputeRadFlux() override;
  void FluxUpdatePseudo1D(); // Helper
  void FluxUpdatePseudo2D(); // Helper

  // --- Member variables ---

  // Scalar output helper functions
  /**
   * @brief Computes Problemspecific Scalar QOI
   *
   */

  double GetCurrentOutflow() override final;
  /**
   * @brief Computes Problemspecific Scalar QOI
   */
  double GetTotalOutflow(unsigned iteration) override final;
  /**
   * @brief Computes Problemspecific Scalar QOI
   */
  double GetMaxOutflow() override final;
  /**
   * @brief Computes Problemspecific Scalar QOI
   */
  double GetFinalTimeAbsorption() override final;
  /**
   * @brief Computes Problemspecific Scalar QOI
   */
  double GetTotalAbsorption() override final;
  /**
   * @brief Computes Problemspecific Scalar QOI
   */
  double GetMaxAbsorption() override final;
  /**
   * @brief Computes Problemspecific Scalar QOI
   */
  double GetTotalAbsorptionCenter() override final;
  /**
   * @brief Computes Problemspecific Scalar QOI
   */
  double GetTotalAbsorptionVertical() override final;
  /**
   * @brief Computes Problemspecific Scalar QOI
   */
  double GetTotalAbsorptionHorizontal() override final;
};

#endif // SNSOLVER_H
