#ifndef SNSOLVER_H
#define SNSOLVER_H

#include "solvers/solverbase.h"

class SNSolver : public Solver
{
  protected:
    Matrix _scatteringKernel; /*!  @brief scattering kernel for the quadrature */

    // quadrature related numbers

    VectorVector _quadPoints; /*!  @brief quadrature points, dim(_quadPoints) = (_nq,spatialDim) */
    Vector _weights;          /*!  @brief quadrature weights, dim(_weights) = (_nq) */

  public:
    /*! @brief SNSolver constructor
     * @param settings stores all needed information
     */
    SNSolver( Config* settings );

  protected:
    // IO
    void virtual PrepareVolumeOutput() override;
    void virtual WriteVolumeOutput( unsigned idx_pseudoTime ) override;

    // Solver
    void virtual FVMUpdate( unsigned idx_energy ) override;
    void virtual FluxUpdate() override;
    void virtual IterPreprocessing( unsigned idx_pseudotime ) override;
    void virtual IterPostprocessing() override;

    // Helper
    void ComputeRadFlux();

    // --- Member variables ---
};

#endif    // SNSOLVER_H
