#ifndef PNSOLVER_H
#define PNSOLVER_H

#include <cmath>

#include "solver.h"
#include "typedef.h"

#include "settings/config.h"

class PNSolver : public Solver
{
  public:
    /**
     * @brief PNSolver constructor
     * @param settings stores all needed information
     */
    PNSolver( Config* settings );
    /**
     * @brief Solve functions runs main time loop
     */
    void Solve() override;
    /**
     * @brief Output solution to VTK file
     */
    void Save() const override;

    void Save( int currEnergy ) const;

  protected:
    // moment orders for P_N
    unsigned _nTotalEntries;    // total number of equations in the system

    VectorVector _sigmaA;    // absorbtion coefficient for all energies
    // System Matrix for x, y and z flux
    //    ==> not needed after computation of A+ and A- ==> maybe safe only temporarly and remove as member?
    SymMatrix _Ax;
    SymMatrix _Ay;
    SymMatrix _Az;

    // Upwinding Matrices
    Matrix _AxPlus;
    Matrix _AxMinus;
    Matrix _AxAbs;
    Matrix _AyPlus;
    Matrix _AyMinus;
    Matrix _AyAbs;
    Matrix _AzPlus;
    Matrix _AzMinus;
    Matrix _AzAbs;

    Vector _scatterMatDiag;    // diagonal of the scattering matrix (its a diagonal matrix by construction)

    // parameter functions for setting up system matrix
    double AParam( int l, int k ) const;
    double BParam( int l, int k ) const;
    double CParam( int l, int k ) const;
    double DParam( int l, int k ) const;
    double EParam( int l, int k ) const;
    double FParam( int l, int k ) const;

    double CTilde( int l, int k ) const;
    double DTilde( int l, int k ) const;
    double ETilde( int l, int k ) const;
    double FTilde( int l, int k ) const;

    // mathematical + index functions
    int Sgn( int k ) const;
    int kPlus( int k ) const;
    int kMinus( int k ) const;

    /*! @brief : computes the global index of the moment corresponding to basis function (l,k)
     *  @param : degree l, it must hold: 0 <= l <=_nq
     *  @param : order k, it must hold: -l <=k <= l
     *  @returns : global index
     */
    int GlobalIndex( int l, int k ) const;

    /*! @brief : Checks, if index invariant for global index holds
     *  @returns : True, if invariant holds, else false
     */
    bool CheckIndex( int l, int k ) const;

    // function for computing and setting up system matrices
    void ComputeSystemMatrices();
    // function for computing and setting up flux matrices for upwinding
    void ComputeFluxComponents();
    // fucntion for computing and setting up diagonal matrix for scattering kernel
    void ComputeScatterMatrix();
    // Computes Legedre polinomial of oder l at point x
    double Legendre( double x, int l );
};

#endif    // PNSOLVER_H
