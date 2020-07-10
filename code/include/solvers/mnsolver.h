#ifndef MNSOLVER_H
#define MNSOLVER_H

#include "solverbase.h"

class EntropyBase;
class QuadratureBase;
class SphericalHarmonics;
class OptimizerBase;

class MNSolver : public Solver
{
  public:
    /**
     * @brief MNSolver constructor
     * @param settings stores all needed information
     */
    MNSolver( Config* settings );

    /*! @brief MNSolver destructor */
    ~MNSolver();

    /**
     * @brief Solve functions runs main time loop
     */
    void Solve() override;
    void Save() const override;                 /*! @brief Save Output solution to VTK file */
    void Save( int currEnergy ) const override; /*! @brief Save Output solution at given energy (pseudo time) to VTK file */

  private:
    unsigned _nTotalEntries;          /*! @brief: Total number of equations in the system */
    unsigned short _nMaxMomentsOrder; /*! @brief: Max Order of Moments */

    VectorVector _sigmaA; /*! @brief: Absorbtion coefficient for all energies */

    SphericalHarmonics* _basis; /*! @brief: Class to compute and store current spherical harmonics basis */
    VectorVector _moments;      /*! @brief: Moment Vector pre-computed at each quadrature point: dim= _nq x _nTotalEntries */

    Vector _scatterMatDiag; /*! @brief: Diagonal of the scattering matrix (its a diagonal matrix by construction) */

    QuadratureBase* _quadrature; /*! @brief: Quadrature rule to compute flux Jacobian */    // Shift this to SolverBase!

    EntropyBase* _entropy; /*! @brief: Class to handle entropy functionals */
    VectorVector _alpha;   /*! @brief: Lagrange Multipliers for Minimal Entropy problem for each gridCell
                                       Layout: _nCells x _nTotalEntries*/

    OptimizerBase* _optimizer; /*! @brief: Class to solve minimal entropy problem */
    /*! @brief : computes the global index of the moment corresponding to basis function (l,k)
     *  @param : degree l, it must hold: 0 <= l <=_nq
     *  @param : order k, it must hold: -l <=k <= l
     *  @returns : global index
     */
    int GlobalIndex( int l, int k ) const;

    /*! @brief : Construct flux by computing the Moment of the  sum of FVM discretization at the interface of cell
     *  @param : idx_cell = current cell id
     *  @returns : sum over all neighbors of flux for all moments at interface of idx_cell, idx_neighbor */
    Vector ConstructFlux( unsigned idx_cell );

    /*! @brief : Pre-Compute Moments at all quadrature points. */
    void ComputeMoments();
};
#endif    // MNSOLVER_H
