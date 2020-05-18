#ifndef NEWTONOPTIMIZER_H
#define NEWTONOPTIMIZER_H

#include <vector>
#include "settings/config.h"

class NewtonOptimizer
{
public:
  NewtonOptimizer();

  /*! \brief: Solves the optimization problem for given moments*/
  Solve(std::vector<std::vector<double>> moments);

  /*! \brief Get pointer to solution vector for each cell */
  std::vector<std::vector<double>> * GetSolution(){return &_lambda};

protected:
  std::vector<std::vector<double>> _lambda; /*! \brief: Solution of the optimization problem */

};

#endif // NEWTONOPTIMIZER_H
