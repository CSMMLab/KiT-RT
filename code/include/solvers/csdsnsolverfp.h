#ifndef CSDSNSOLVERFP_H
#define CSDSNSOLVERFP_H

#include "icru.h"
#include "solvers/snsolver.h"

class Physics;

class CSDSNSolverFP : public SNSolver
{
  private:
    std::vector<double> _dose; /*! @brief: TODO */

    // Physics acess
    Vector _energies; /*! @brief: energy levels for CSD, lenght = _nEnergies */
    Vector _angle;    /*! @brief: angles for SN */
    std::vector<double> _density;  /*! @brief: patient density for each grid cell */

    std::vector<Matrix> _sigmaSE; /*!  @brief scattering cross section for all energies*/
    Vector _sigmaTE;              /*!  @brief total cross section for all energies*/

    Matrix _L; /*!  @brief Laplace Beltrami Matrix */
    Matrix _IL; /*!  @brief Laplace Beltrami Matrix */

    double _alpha;
    double _beta;

    Vector _xi1;
    Vector _xi2;
    Matrix _xi;

    bool _RT;

    double _energyMin;
    double _energyMax;

  public:
    /**
     * @brief CSDSNSolverFP constructor
     * @param settings stores all needed information
     */
    CSDSNSolverFP( Config* settings );
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

#endif    // CSDSNSOLVERFP_H
