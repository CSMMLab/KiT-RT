#ifndef SNSOLVER_H
#define SNSOLVER_H

#include "solvers/solverbase.h"

class SNSolver : public Solver
{
  protected:
    Matrix _scatteringKernel; /*!  @brief scattering kernel for the quadrature */

    // quadrature related numbers

    VectorVector _quadPoints; /*!  @brief quadrature points, dim(_quadPoints) = (_nq,spatialDim) */
    Vector _weights;          /*!  @brief quadrature weights, dim(_weights) = (_nq) */

  public:
    /**
     * @brief SNSolver constructor
     * @param settings stores all needed information
     */
    SNSolver( Config* settings );
    /**
     * @brief Solve functions runs main time loop
     */
    void Solve() override;
    /**
     * @brief Output solution to VTK file
     */
    void Save() const override;
    void Save( int currEnergy ) const override;
};

#endif    // SNSOLVER_H
