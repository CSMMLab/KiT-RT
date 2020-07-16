#ifndef CSDSNSOLVER_H
#define CSDSNSOLVER_H

#include "solvers/snsolver.h"

class CSDSNSolver : public SNSolver
{
  private:
    std::vector<double> _dose;
    std::vector<double> _density;    // patient density, dim(_density) = _nCells

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
