/*!
 * \file GlobalConstants.h
 * \brief All global defined (physical) constants, enums etc
 * \author <blank>
 * \version 0.0
 *
 */

#ifndef GLOBAL_CONSTANTS_H
#define GLOBAL_CONSTANTS_H

#include <cmath>
#include <map>
#include <string>

// --- Definition for global constants goes here ---

const double PI_NUMBER             = M_PI; /*!< @brief Pi number. */
const unsigned int MAX_STRING_SIZE = 200;  /*!< @brief Maximum size for strings. */

// --- Definition of enums goes here ---

enum BOUNDARY_TYPE { DIRICHLET, NEUMANN, NONE, INVALID };

// --- Definition of enums for EnumOptions goes here ---

/*! @brief Enum for all currently available quadratures in rtsn.
 *         Option enums are written in capital letters with underscores as spaces (e.g option "time integration" has option enum "TIME_INTEGRATION")
 */
enum QUAD_NAME { QUAD_MonteCarlo, QUAD_GaussLegendreTensorized, QUAD_GaussLegendre1D, QUAD_LevelSymmetric, QUAD_Lebedev, QUAD_LDFESA };

/*! @brief Conversion Map String to enum
 */
inline std::map<std::string, QUAD_NAME> Quadrature_Map{ { "MONTE_CARLO", QUAD_MonteCarlo },
                                                        { "GAUSS_LEGENDRE_TENSORIZED", QUAD_GaussLegendreTensorized },
                                                        { "GAUSS_LEGENDRE_1D", QUAD_GaussLegendre1D },
                                                        { "LEVEL_SYMMETRIC", QUAD_LevelSymmetric },
                                                        { "LEBEDEV", QUAD_Lebedev },
                                                        { "LDFESA", QUAD_LDFESA } };

// Problem name
enum PROBLEM_NAME { PROBLEM_LineSource, PROBLEM_Checkerboard, PROBLEM_ElectronRT, PROBLEM_WaterPhantom, PROBLEM_LineSource_Pseudo_1D };

inline std::map<std::string, PROBLEM_NAME> Problem_Map{ { "LINESOURCE", PROBLEM_LineSource },
                                                        { "CHECKERBOARD", PROBLEM_Checkerboard },
                                                        { "ELECTRONRT", PROBLEM_ElectronRT },
                                                        { "WATERPHANTOM", PROBLEM_WaterPhantom },
                                                        { "LINESOURCE_PSEUDO_1D", PROBLEM_LineSource_Pseudo_1D } };

// Kernel name
enum KERNEL_NAME { KERNEL_Isotropic, KERNEL_Isotropic1D };

inline std::map<std::string, KERNEL_NAME> Kernel_Map{ { "ISOTROPIC", KERNEL_Isotropic }, { "ISOTROPIC_1D", KERNEL_Isotropic1D } };

// Solver name
enum SOLVER_NAME { SN_SOLVER, PN_SOLVER, MN_SOLVER };

inline std::map<std::string, SOLVER_NAME> Solver_Map{ { "SN_SOLVER", SN_SOLVER }, { "PN_SOLVER", PN_SOLVER }, { "MN_SOLVER", MN_SOLVER } };

// Entropy functional
enum ENTROPY_NAME { QUADRATIC, MAXWELL_BOLZMANN, BOSE_EINSTEIN, FERMI_DIRAC };

inline std::map<std::string, ENTROPY_NAME> Entropy_Map{
    { "QUADRATIC", QUADRATIC }, { "MAXWELL_BOLZMANN", MAXWELL_BOLZMANN }, { "BOSE_EINSTEIN", BOSE_EINSTEIN }, { "FERMI_DIRAC", FERMI_DIRAC } };

// Otpimizer
enum OPTIMIZER_NAME { NEWTON, ML };

inline std::map<std::string, OPTIMIZER_NAME> Optimizer_Map{ { "NEWTON", NEWTON }, { "ML", ML } };
#endif    // GLOBAL_CONSTANTS_H
