#ifndef MNSOLVER_H
#define MNSOLVER_H

#include "solverbase.h"
#include "sphericalharmonics.h"

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

    VectorVector _sigmaA;        /*! @brief: Absorbtion coefficient for all energies */
    SphericalHarmonics _basis;   /*! @brief: Class to compute and store current spherical harmonics basis */
    Vector _scatterMatDiag;      /*! @brief: Diagonal of the scattering matrix (its a diagonal matrix by construction) */
    VectorVector _A;             /*! @brief: Diagonals of the system matrices  (each sysmatric is a diagonal matrix by construction)
                                             layout: Nx2 (in 2d), N = size of system */
    QuadratureBase* _quadrature; /*! @brief: Quadrature rule to compute flux Jacobian */

    /*! @brief : computes the global index of the moment corresponding to basis function (l,k)
     *  @param : degree l, it must hold: 0 <= l <=_nq
     *  @param : order k, it must hold: -l <=k <= l
     *  @returns : global index
     */
    int GlobalIndex( int l, int k ) const;

    /*! @brief : function for computing and setting up the System Matrices.
                 The System Matrices are diagonal, so there is no need to diagonalize*/
    void ComputeSystemMatrices();
};
#endif    // MNSOLVER_H
