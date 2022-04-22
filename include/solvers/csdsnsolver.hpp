#ifndef CSDSNSOLVER_H
#define CSDSNSOLVER_H

// externals
#include "spdlog/spdlog.h"
#include <mpi.h>

#include "common/config.hpp"
#include "common/io.hpp"
#include "fluxes/numericalflux.hpp"
#include "kernels/scatteringkernelbase.hpp"
#include "problems/problembase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "solvers/snsolver.hpp"
#include "velocitybasis/sphericalharmonics.hpp"

class Physics;

class CSDSNSolver : public SNSolver
{
  private:
    std::vector<double> _dose; /*!< @brief TODO */

    Matrix _O;    // Transformation Matrix from Moments to Ordinates
    Matrix _M;    // Transformation Matrix from Ordinates to Moments

    VectorVector _quadPointsSphere;
    Vector _eTrafo;

    // Helper variables
    unsigned _polyDegreeBasis;

  public:
    /**
     * @brief CSDSolverTrafoFP2D constructor
     * @param settings stores all needed information
     */
    CSDSNSolver( Config* settings );

    virtual ~CSDSNSolver() {}

    /**
     * @brief Output solution to VTK file (LEGACY. But not remove, to compare)
     */
    // virtual void Save() const;
    // virtual void Save( int currEnergy ) const;

  private:
    void GenerateEnergyGrid( bool refinement );

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

#endif    // CSDSNSOLVER_H
