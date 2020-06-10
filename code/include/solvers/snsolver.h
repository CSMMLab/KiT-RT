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
    virtual void Solve();
    /**
     * @brief Output solution to VTK file
     */
    virtual void Save() const;
    virtual void Save( int currEnergy ) const;
};

#endif    // SNSOLVER_H
