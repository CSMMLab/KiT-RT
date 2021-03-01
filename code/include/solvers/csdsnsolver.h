#ifndef CSDSNSOLVER_H
#define CSDSNSOLVER_H

#include "common/typedef.h"    // for Matrix, Vector
#include "solvers/snsolver.h"

#include <vector>    // for vector

class Config;

class CSDSNSolver : public SNSolver
{
  private:
    std::vector<double> _dose; /*!< @brief TODO */

    // Physics acess
    Vector _energies; /*!< @brief energy levels for CSD, lenght = _nEnergies */
    Vector _angle;    /*!< @brief angles for SN */

    std::vector<Matrix> _sigmaSE; /*!<  @brief scattering cross section for all energies*/
    Vector _sigmaTE;              /*!<  @brief total cross section for all energies*/

  public:
    /**
     * @brief CSDSNSolver constructor
     * @param settings stores all needed information
     */
    CSDSNSolver( Config* settings );

    virtual ~CSDSNSolver() {}

    /**
     * @brief Solve functions runs main time loop
     */
    void Solve() override;

  private:
    void PrepareVolumeOutput() override;
    void WriteVolumeOutput( unsigned idx_pseudoTime ) override;
};

#endif    // CSDSNSOLVER_H
