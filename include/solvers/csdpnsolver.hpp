#ifndef CSDPNSOLVER_H
#define CSDPNSOLVER_H

#include "solvers/pnsolver.hpp"

class SphericalBase;

class CSDPNSolver : public PNSolver
{
  public:
    /**
     * @brief CSDPNSolver constructor
     * @param settings stores all needed information
     */
    CSDPNSolver( Config* settings );

    virtual ~CSDPNSolver();

  private:
    std::vector<double> _dose; /*!< @brief Radiation Dose */
    Vector _eTrafo;            /*!< @brief Transformed energy grid */
    Vector _sigmaTAtEnergy;    /*!< @brief Scattercoefficient at energy grid */

    void SolverPreprocessing() override;
    void IterPreprocessing( unsigned idx_iter ) override;
    void IterPostprocessing( unsigned idx_iter ) override;
    void FluxUpdate() override;
    void FVMUpdate( unsigned idx_energy ) override;
    void PrepareVolumeOutput() override;
    void WriteVolumeOutput( unsigned idx_pseudoTime ) override;

    // Helper Functions for CSD
    // double NormPDF( double x, double mu, double sigma );
};

#endif    // CSDPNSOLVER_H
