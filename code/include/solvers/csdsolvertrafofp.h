#ifndef CSDSOLVERTRAFOFP_H
#define CSDSOLVERTRAFOFP_H

#include "icru.h"
#include "solvers/snsolver.h"

class Physics;

class CSDSolverTrafoFP : public SNSolver
{
  private:
    std::vector<double> _dose; /*! @brief: TODO */

    // Physics acess
    Vector _energies; /*! @brief: energy levels for CSD, lenght = _nEnergies */
    Vector _density;  /*! @brief: patient density for each grid cell */

    std::vector<Matrix> _sigmaSE; /*!  @brief scattering cross section for all energies*/
    Vector _sigmaTE;              /*!  @brief total cross section for all energies*/

    Matrix _L;  /*!  @brief Laplace Beltrami Matrix */
    Matrix _IL; /*!  @brief Laplace Beltrami Matrix */

    double _alpha;
    double _alpha2;
    double _beta;

    Vector _xi1;
    Vector _xi2;
    Matrix _xi;

    bool _RT;

    double _energyMin;
    double _energyMax;

    void GenerateEnergyGrid( bool refinement );

  public:
    /**
     * @brief CSDSolverTrafoFP constructor
     * @param settings stores all needed information
     */
    CSDSolverTrafoFP( Config* settings );
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

#endif    // CSDSOLVERTRAFOFP_H
