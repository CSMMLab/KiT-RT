#ifndef CSDPNSOLVER_H
#define CSDPNSOLVER_H

#include "solvers/pnsolver.h"

class CSDPNSolver : public PNSolver
{
  private:
    std::vector<double> _dose; /*!< @brief TODO */

    // Physics acess
    Vector _energies; /*!< @brief energy levels for CSD, lenght = _nEnergies */
    Vector _angle;    /*!< @brief angles */

    std::vector<Matrix> _sigmaSE; /*!<  @brief scattering cross section for all energies*/
    Vector _sigmaTE;              /*!<  @brief total cross section for all energies*/

  public:
    /**
     * @brief CSDPNSolver constructor
     * @param settings stores all needed information
     */
    CSDPNSolver( Config* settings );

    virtual ~CSDPNSolver() {}

    /**
     * @brief Solve functions runs main time loop
     */
    void Solve() override;

  private:
    void PrepareVolumeOutput() override;
    void WriteVolumeOutput( unsigned idx_pseudoTime ) override;
};

#endif    // CSDPNSOLVER_H
