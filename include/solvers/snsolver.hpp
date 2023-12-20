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
    SNSolver( Config* settings );

    virtual ~SNSolver() {}

  protected:
    // IO
    void virtual PrepareVolumeOutput() override;
    void virtual WriteVolumeOutput( unsigned idx_iter ) override;

    // Solver
    void virtual FVMUpdate( unsigned idx_iter ) override;
    void virtual FluxUpdate() override;
    void virtual IterPreprocessing( unsigned idx_iter ) override;

    // Helper
    void ComputeScalarFlux() override;
    void FluxUpdatePseudo1D();    // Helper
    void FluxUpdatePseudo2D();    // Helper
};

#endif    // SNSOLVER_H
