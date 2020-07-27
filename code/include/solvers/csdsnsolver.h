#ifndef CSDSNSOLVER_H
#define CSDSNSOLVER_H

#include "solvers/snsolver.h"

class Physics;

class CSDSNSolver : public SNSolver
{
  private:
    std::vector<double> _dose; /*! @brief: TODO */

    // Physics acess
    Vector _energies; /*! @brief: energy levels for CSD, lenght = _nEnergies */
    Vector _angle;    /*! @brief: angles for SN */
    Vector _density;  /*! @brief: patient density for each grid cell */

    VectorVector _sigmaSE; /*!  @brief scattering cross section for all energies*/
    Vector _sigmaTE;       /*!  @brief total cross section for all energies*/

  public:
    /**
     * @brief CSDSNSolver constructor
     * @param settings stores all needed information
     */
    CSDSNSolver( Config* settings );
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

#endif    // CSDSNSOLVER_H
