#ifndef CSDSOLVERTRAFOFP_H
#define CSDSOLVERTRAFOFP_H

#include "solvers/snsolver.hpp"

class Physics;

class CSDSolverTrafoFP : public SNSolver
{
  private:
    std::vector<double> _dose; /*!< @brief TODO */

    // Helper Variables

    Matrix _L;  /*!<  @brief Laplace Beltrami Matrix */
    Matrix _IL; /*!<  @brief Laplace Beltrami Matrix */

    double _alpha;  /*!<  @brief  Coefficient of GFP operators (see Olbrant 2010, Appendix B)*/
    double _alpha2; /*!<  @brief  Coefficient of GFP operators (see Olbrant 2010, Appendix B)*/
    double _beta;   /*!<  @brief  Coefficient of GFP operators (see Olbrant 2010, Appendix B)*/

    Matrix _xi; /*!<  @brief matrix of transport coefficients */
    Vector _xi1;
    Vector _xi2;

    unsigned _FPMethod; /*!<  @brief Encodes different ways of computing coefficients alpha, alpha2 & beta, _FPMethod == 1, 2 ,3 stand for methods
                           with increasing accuracy (see Olbrant 2010, Appendix B)*/

    bool _RT; /*!<  @brief radiotherapy application (on/off), if true use crosssections + stopping powers from database  */

    double _energyMin; /*!<  @brief minimal energy in energy grid*/
    double _energyMax; /*!<  @brief maximal energy in energy grid*/

    void GenerateEnergyGrid( bool refinement );

    // Helper variables
    Vector _energiesOrig; /*!< @brief original energy levels for CSD, lenght = _nEnergies */
    Matrix _identity;     /*!< @brief: identity matrix for FP scattering. Dim (_nq,_nq)*/

  public:
    /**
     * @brief CSDSolverTrafoFP constructor
     * @param settings stores all needed information
     */
    CSDSolverTrafoFP( Config* settings );

    virtual ~CSDSolverTrafoFP() {}

  private:
    // IO
    void PrepareVolumeOutput() override final;
    void WriteVolumeOutput( unsigned idx_pseudoTime ) override final;

    // Solver
    void FVMUpdate( unsigned idx_energy ) override final;
    void FluxUpdate() override final;
    void IterPreprocessing( unsigned idx_pseudotime ) override final;
    void virtual IterPostprocessing( unsigned idx_pseudotime ) override final;
    void SolverPreprocessing() override final;
};

#endif    // CSDSOLVERTRAFOFP_H
