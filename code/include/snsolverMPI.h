#ifndef SNSOLVERMPI_H
#define SNSOLVERMPI_H

#include <mpi.h>

#include "solver.h"
#include "typedef.h"

#include "settings/config.h"

class SNSolverMPI : public Solver
{
  public:
    /**
     * @brief SNSolver constructor
     * @param settings stores all needed information
     */
    SNSolverMPI( Config* settings );
    /**
     * @brief Solve functions runs main time loop
     */
    virtual void Solve();
    /**
     * @brief Output solution to VTK file
     */
    virtual void Save() const;
};

#endif    // SNSOLVERMPI_H
