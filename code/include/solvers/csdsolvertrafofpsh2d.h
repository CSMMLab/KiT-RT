#ifndef CSDSOLVERTRAFOFPSH2D_H
#define CSDSOLVERTRAFOFPSH2D_H

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

#include "common/config.h"
#include "common/io.h"
#include "fluxes/numericalflux.h"
#include "kernels/scatteringkernelbase.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"
#include "solvers/snsolver.h"
#include "toolboxes/sphericalharmonics.h"

class Physics;

class CSDSolverTrafoFPSH2D : public SNSolver
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

    Matrix _O;
    Matrix _S;
    Matrix _M;

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
    CSDSolverTrafoFPSH2D( Config* settings );
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

#endif    // CSDSOLVERTRAFOFPSH2D_H
