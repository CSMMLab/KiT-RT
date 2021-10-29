#ifndef CSDSNSOLVERFP_H
#define CSDSNSOLVERFP_H

#include "solvers/snsolver.hpp"

class Physics;

class CSDSNSolverFP : public SNSolver
{
  private:
    std::vector<double> _dose; /*!< @brief TODO */

    // Physics acess
    Vector _angle; /*!< @brief angles for SN */

    Matrix _L;  /*!<  @brief Laplace Beltrami Matrix */
    Matrix _IL; /*!<  @brief Laplace Beltrami Matrix */

    double _alpha; /*!<  @brief  Coefficient of GFP operators (see Olbrant 2010, eq. (8)*/
    double _beta;  /*!<  @brief  Coefficient of GFP operators (see Olbrant 2010, eq. (8)*/

    Matrix _xi; /*!<  @brief matrix of transport coefficients */

    bool _RT; /*!<  @brief radiotherapy application (on/off), if true use crosssections + stopping powers from database  */

    double _energyMin; /*!<  @brief minimal energy in energy grid*/
    double _energyMax; /*!<  @brief maximal energy in energy grid*/

  public:
    /**
     * @brief CSDSNSolverFP constructor
     * @param settings stores all needed information
     */
    CSDSNSolverFP( Config* settings );

    virtual ~CSDSNSolverFP() {}

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
