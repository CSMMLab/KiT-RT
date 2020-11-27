#ifndef SNSOLVERMPI_H
#define SNSOLVERMPI_H

#include <mpi.h>

#include "solvers/snsolver.h"

class SNSolverMPI : public SNSolver
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
    virtual void PrintVolumeOutput() const;
};

#endif    // SNSOLVERMPI_H
