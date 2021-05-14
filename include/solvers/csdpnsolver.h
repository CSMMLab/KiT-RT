#ifndef CSDPNSOLVER_H
#define CSDPNSOLVER_H

#include "solvers/pnsolver.h"

class CSDPNSolver : public PNSolver
{
  private:
    std::vector<double> _dose; /*!< @brief TODO */

    // Physics acess
    Vector _angle; /*!< @brief angles */

    std::vector<Matrix> _sigmaSE; /*!<  @brief scattering cross section for all energies*/
    Vector _sigmaTE;              /*!<  @brief total cross section for all energies*/

  public:
    /**
     * @brief CSDPNSolver constructor
     * @param settings stores all needed information
     */
    CSDPNSolver( Config* settings );

    virtual ~CSDPNSolver() {}

    // virtual Solve() override;

  private:
    void SolverPreprocessing() override;

    void IterPreprocessing( unsigned /*idx_iter*/ ) override;
    void IterPostprocessing( unsigned /*idx_iter*/ ) override;

    void FluxUpdate() override;
    void FVMUpdate( unsigned idx_energy ) override;

    void PrepareVolumeOutput() override;
    void WriteVolumeOutput( unsigned idx_pseudoTime ) override;
};

#endif    // CSDPNSOLVER_H
