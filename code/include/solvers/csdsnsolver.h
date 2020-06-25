#ifndef CSDSNSOLVER_H
#define CSDSNSOLVER_H

#include <mpi.h>

#include "solvers/solverbase.h"

class CSDSNSolver : public Solver
{
  private:
    std::vector<double> _dose;

  public:
    /**
     * @brief CSDSNSolver constructor
     * @param settings stores all needed information
     */
    CSDSNSolver( Config* settings );
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

#endif    // CSDSNSOLVER_H
