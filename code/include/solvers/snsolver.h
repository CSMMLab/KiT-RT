#ifndef SNSOLVER_H
#define SNSOLVER_H

#include <mpi.h>

#include "solvers/solverbase.h"

class SNSolver : public Solver
{
  private:
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
