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
    Vector _angle;    /*! @brief: angles for SN */

    std::vector<Matrix> _sigmaSE; /*!  @brief scattering cross section for all energies*/
    Vector _sigmaTE;              /*!  @brief total cross section for all energies*/

    Matrix _L;  /*!  @brief Laplace Beltrami Matrix */
    Matrix _IL; /*!  @brief Laplace Beltrami Matrix */

    double _alpha; /*!  @brief  Coefficient of GFP operators (see Olbrant 2010, Appendix B)*/
    double _alpha2; /*!  @brief  Coefficient of GFP operators (see Olbrant 2010, Appendix B)*/
    double _beta; /*!  @brief  Coefficient of GFP operators (see Olbrant 2010, Appendix B)*/

    Matrix _xi;  /*!  @brief matrix of transport coefficients */
    Vector _xi1;
    Vector _xi2;
    
    unsigned _FPMethod; /*!  @brief Encodes different ways of computing coefficients alpha, alpha2 & beta, _FPMethod == 1, 2 ,3 stand for methods with increasing accuracy (see Olbrant 2010, Appendix B)*/

    bool _RT; /*!  @brief radiotherapy application (on/off), if true use crosssections + stopping powers from database  */

    double _energyMin; /*!  @brief minimal energy in energy grid*/
    double _energyMax; /*!  @brief maximal energy in energy grid*/

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