#ifndef CSDSOLVERTRAFOFP2D_H
#define CSDSOLVERTRAFOFP2D_H

#include "solvers/snsolver.hpp"

class Physics;

class CSDSolverTrafoFP2D : public SNSolver
{
  private:
    std::vector<double> _dose; /*!< @brief TODO */
    Vector _sigmaTAtEnergy;
    unsigned short _polyDegreeBasis; /*!< @brief Max polynomial degree of the basis */

    // Helper Variables

    Matrix _L;  /*!<  @brief Laplace Beltrami Matrix */
    Matrix _IL; /*!<  @brief Laplace Beltrami Matrix */

    VectorVector _quadPointsSphere;
    Vector _mu;
    Vector _phi;
    Vector _wp;
    Vector _wa;

    double _alpha;
    double _alpha2;
    double _beta;

    double _E_cutoff;
    Vector _eTrafo;

    Vector _xi1;
    Vector _xi2;
    Matrix _xi;

    unsigned _FPMethod;

    bool _RT;

    double _energyMin;
    double _energyMax;

    // Helper variables
    Vector _energiesOrig; /*!< @brief original energy levels for CSD, lenght = _nEnergies */
    Matrix _identity;     /*!< @brief: identity matrix for FP scattering. Dim (_nq,_nq)*/
    double _densityMin;   /*!< @brief Minimal density of _density vector */

  public:
    /**
     * @brief CSDSolverTrafoFP2D constructor
     * @param settings stores all needed information
     */
    CSDSolverTrafoFP2D( Config* settings );

    virtual ~CSDSolverTrafoFP2D() {}

    /**
     * @brief Output solution to VTK file
     */

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
    double NormPDF( double x, double mu, double sigma );
};

#endif    // CSDSOLVERTRAFOFP2D_H
