#ifndef CSDSOLVERTRAFOFP2D_H
#define CSDSOLVERTRAFOFP2D_H

#include "solvers/snsolver.h"

class Physics;

class CSDSolverTrafoFP2D : public SNSolver
{
  private:
    std::vector<double> _dose; /*! @brief: TODO */

    // Physics acess
    Vector _energies; /*! @brief: energy levels for CSD, lenght = _nEnergies */
    Vector _angle;    /*! @brief: angles for SN */

    std::vector<Matrix> _sigmaSE; /*!  @brief scattering cross section for all energies*/
    Vector _sigmaTE;              /*!  @brief total cross section for all energies*/

    Matrix _L;  /*!  @brief Laplace Beltrami Matrix */
    Matrix _IL; /*!  @brief Laplace Beltrami Matrix */

    VectorVector _quadPoints;
    VectorVector _quadPointsSphere;
    Vector _weights;
    Vector _mu;
    Vector _phi;
    Vector _wp;
    Vector _wa;

    double _alpha;
    double _alpha2;
    double _beta;

    Vector _xi1;
    Vector _xi2;
    Matrix _xi;

    unsigned _FPMethod;

    bool _RT;

    double _energyMin;
    double _energyMax;

    void GenerateEnergyGrid( bool refinement );

  public:
    /**
     * @brief CSDSolverTrafoFP2D constructor
     * @param settings stores all needed information
     */
    CSDSolverTrafoFP2D( Config* settings );
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

#endif    // CSDSOLVERTRAFOFP2D_H
