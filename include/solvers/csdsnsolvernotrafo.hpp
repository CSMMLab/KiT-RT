#ifndef CSDSNSOLVERNOTRAFO_H
#define CSDSNSOLVERNOTRAFO_H

#include "solvers/snsolver.hpp"

class Physics;

class CSDSNSolverNoTrafo : public SNSolver
{
  private:
    std::vector<double> _dose; /*!< @brief TODO */

    // Physics acess
    Vector _angle; /*!< @brief angles for SN */

    std::vector<Matrix> _sigmaSE; /*!<  @brief scattering cross section for all energies*/
    Vector _sigmaTE;              /*!<  @brief total cross section for all energies*/

  public:
    /**
     * @brief CSDSNSolverNoTrafo constructor
     * @param settings stores all needed information
     */
    CSDSNSolverNoTrafo( Config* settings );

    virtual ~CSDSNSolverNoTrafo() {}

    /**
     * @brief Solve functions runs main time loop
     */
    virtual void Solve();
    /**
     * @brief Output solution to VTK file
     */
    virtual void Save() const;
    virtual void Save( int currEnergy ) const;
};

#endif    // CSDSNSOLVERNOTRAFO_H
