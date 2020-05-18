/*!
 * \file advectionsolver.h
 * \brief solver for entropy closure of advection equation.
 * \author S. Schotthöfer
 */

#ifndef ADVECTIONSOLVER_H
#define ADVECTIONSOLVER_H

#include <vector>
#include "settings/config.h"


class AdvectionSolver{
public:
  AdvectionSolver(Config* settings);

  /*! \brief: Solves the advection problem for given lambda */
  void Solve(std::vector<std::vector<double>> lambda);

  /*! \brief: Reconstructs kintec density from given lambda and legendre dual of the entropy functional */
  std::vector<double> * ReconstructKineticDensity(std::vector<std::vector<double>> lambda);

  /* ---- Getter ---- */
  /*! \brief Get pointer to solution vector for each cell (moments) */
  std::vector<std::vector<double>> * GetSolution(void){return &_moments};


protected:
  std::vector<std::vector<double>> _moments; /*! \brief: Solution of the moment advection problem */
  std::vector<std::vector<double>> _kineticDensity; /*! \brief: Reconstructed kinetic Density from given lamnda */

}

#endif // ADVECTIONSOLVER_H_H
