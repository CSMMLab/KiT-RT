#ifndef MNSOLVER_H
#define MNSOLVER_H

#include <pnsolver.h>
#include <sphericalharmonics.h>

class MNSolver : Solver
{
  public:
    /**
     * @brief MNSolver constructor
     * @param settings stores all needed information
     */
    MNSolver( Config* settings );

    /**
     * @brief Solve functions runs main time loop
     */
    void Solve() override;

  private:
    /*! @brief: Total number of equations in the system */
    unsigned _nTotalEntries;
    /*! @brief: Max Order of Moments */
    unsigned short _nMaxMomentsOrder;

    /*! @brief: Absorbtion coefficient for all energies */
    VectorVector _sigmaA;

    /*! @brief: Class to compute and store current spherical harmonics basis */
    SphericalHarmonics _basis;

    /*! @brief: Diagonal of the scattering matrix (its a diagonal matrix by construction) */
    Vector _scatterMatDiag;

    /*! @brief: Diagonals of the system matrices  (each sysmatric is a diagonal matrix by construction)
                layout: Nx2 (in 2d), N = size of system */
    VectorVector _A;
    ///*! @brief: Diagonal of the system matrix in x direction (its a diagonal matrix by construction) */
    // Vector _Ax;
    //
    ///*! @brief: Diagonal of the system matrix in y direction (its a diagonal matrix by construction) */
    // Vector _Ay;
    //
    ///*! @brief: Diagonal of the system matrix in z direction (its a diagonal matrix by construction) */
    // Vector _Az;

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
