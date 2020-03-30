#ifndef SNSOLVER_H
#define SNSOLVER_H

#include "solver.h"
#include "typedef.h"

class SNSolver : public Solver
{
public:
    /**
     * @brief SNSolver constructor
     * @param settings stores all needed information
     */
    SNSolver(Settings* settings);
    /**
     * @brief Solve functions runs main time loop
     */
    virtual void Solve();
};

#endif // SNSOLVER_H
