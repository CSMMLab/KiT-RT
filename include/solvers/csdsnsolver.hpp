#ifndef CSDSNSOLVER_H
#define CSDSNSOLVER_H

#include "solvers/snsolver.hpp"

class CSDSNSolver : public SNSolver
{
  private:
    std::vector<double> _dose; /*!< @brief TODO */

    // Physics acess
    Vector _angle; /*!< @brief angles for SN */

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
