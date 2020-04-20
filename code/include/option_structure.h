/*!
 * \file option_structure.h
 * \brief All global defined (physical) constants, enums etc
 * \author <blank>
 * \version 0.0
 *
 */

#ifndef OPTION_STRUCTURE_H
#define OPTION_STRUCTURE_H

#include <map>
#include <cmath>

// Definition for global constants goes here
const double PI_NUMBER = 4.0 * atan(1.0);  /*!< \brief Pi number. */


// Definition of enums goes here
/*! \brief Enum for all currently available quadratures in rtsn.
 *         Option enums are written in capital letters with underscores as spaces (e.g option "time integration" has option enum "TIME_INTEGRATION")
 */
enum QUAD_NAME{
    QUAD_MonteCarlo,
    QUAD_GaussLegendreTensorized,
    QUAD_LevelSymmetric,
    QUAD_Lebedev,
    QUAD_LDFESA
};

/*! \brief Conversion Map String to enum
 */
inline std::map <std::string,QUAD_NAME> Quadrature_Map {
  { "MONTE_CARLO", QUAD_MonteCarlo },
  { "GAUSS_LEGENDRE_TENSORIZED", QUAD_GaussLegendreTensorized },
  { "LEVEL_SYMMETRIC", QUAD_LevelSymmetric },
  { "LEBEDEV", QUAD_Lebedev },
  { "LDFESA",QUAD_LDFESA }
};

#endif // OPTION_STRUCTURE_H
