#ifndef CSDMNSolver_H
#define CSDMNSolver_H

#include "solvers/mnsolver.hpp"

class CSDMNSolver : public MNSolver
{
  public:
    /**
     * @brief CSDMNSolver constructor
     * @param settings stores all needed information
     */
    CSDMNSolver( Config* settings );

    virtual ~CSDMNSolver();

  private:
    std::vector<double> _dose; /*!< @brief Radiation Dose */
    Vector _eTrafo;            /*!< @brief Transformed energy grid */
    Vector _sigmaTAtEnergy;    /*!< @brief Scattercoefficient at energy grid */

    void SolverPreprocessing() override;
    void IterPreprocessing( unsigned /*idx_iter*/ ) override;
    void IterPostprocessing( unsigned /*idx_iter*/ ) override;
    void FluxUpdate() override;
    void FVMUpdate( unsigned idx_energy ) override;
    void PrepareVolumeOutput() override;
    void WriteVolumeOutput( unsigned idx_pseudoTime ) override;
    Vector ConstructFlux( unsigned idx_cell ) override;

    // Helper Functions for CSD
    double NormPDF( double x, double mu, double sigma );
    // Vector Time2Energy( const Vector& t, const double E_CutOff );
    // double Time2Energy( const double t, const double E_CutOff );
    // Vector Energy2Time( const Vector& E, const double E_CutOff );
    // double Energy2Time( const double E, const double E_CutOff );
};

#endif    // CSDMNSolver_H
