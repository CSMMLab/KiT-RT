#ifndef SNSOLVER_H
#define SNSOLVER_H

#include "solver.h"
#include "typedef.h"

#include "settings/CConfig.h"

class SNSolver : public Solver
{
  public:
    /**
     * @brief SNSolver constructor
     * @param settings stores all needed information
     */
    SNSolver( CConfig* settings );
    /**
     * @brief Solve functions runs main time loop
     */
    virtual void Solve();
    /**
     * @brief Output solution to VTK file
     */
    virtual void Save() const;
    void SolveMPI();    // can be deleated later
};

#endif    // SNSOLVER_H
