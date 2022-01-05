/*!
 * \file globalconstants.h
 * \brief All global defined (physical) constants, enums etc
 * \version 0.0
 *
 */

#ifndef GLOBAL_CONSTANTS_H
#define GLOBAL_CONSTANTS_H

#include <cmath>
#include <map>
#include <string>
#include <vector>

// --- Definition for global constants goes here ---

static const long double PI           = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
static const long double EULER        = 2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274;
static const long double INV_SQRT_2PI = 1.0 / std::sqrt( 2.0 * PI );

static const long double AVOGADRO_CONSTANT               = 6.02214129e23;       // 1/mol
static const long double ELECTRON_MASS                   = 9.10938291e-31;      // kg
static const long double ELECTRON_ENERGY                 = 0.51099895;          // MeV
static const long double ELECTRON_RADIUS                 = 2.8179403262e-15;    // m
static const long double ELECTRON_MASS_ENERGY_EQUIVALENT = 8.18710506e-14;      // J
static const long double ELECTRON_VOLT                   = 1.602176565e-19;     // J
static const long double ELECTRIC_CONSTANT               = 8.854187817e-12;     // F/m
static const long double ELEMENTARY_CHARGE               = 1.602176565e-19;     // C
static const long double PLANCK_CONSTANT                 = 6.62606957e-34;      // J s
static const long double BOHR_RADIUS                     = 5.2917721092e-11;    // m

static const std::vector<long double> H2OMassFractions{ 0.111894, 0.888106 };
static const std::vector<long double> H2OAtomicNumbers{ 1.0, 8.0 };
static const std::vector<long double> H2OAtomicWeights{ 1.008, 15.999 };
static const long double H2OAtomicWeight     = 2 * H2OAtomicWeights[0] + H2OAtomicWeights[1];
static const long double H2OElectronDensity  = 3.3428847e23;                                            // 1/cm3
static const long double H2OMassDensity      = 1.0;                                                     // g/cm3
static const long double H2OMolecularDensity = AVOGADRO_CONSTANT * H2OMassDensity / H2OAtomicWeight;    // 1/cm3

static const long double C1      = 2.0 * PI * std::pow( ELEMENTARY_CHARGE, 4u ) / std::pow( 4.0 * PI * ELECTRIC_CONSTANT, 2u );
static const long double C2      = C1 / 16.0;
static const long double SQRT05E = std::sqrt( 0.5 * EULER );
static const long double aH = PLANCK_CONSTANT * PLANCK_CONSTANT * ELECTRIC_CONSTANT / ( PI * ELECTRON_MASS * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE );

const unsigned int MAX_STRING_SIZE = 200; /*!< @brief Maximum size for strings. */

// --- Definition of enums goes here ---

enum BOUNDARY_TYPE { DIRICHLET, NEUMANN, NONE, INVALID };

// --- Definition of enums for EnumOptions goes here ---

/*! @brief Enum for all currently available quadratures in rtsn.
 *         Option enums are written in capital letters with underscores as spaces (e.g option "time integration" has option enum "TIME_INTEGRATION")
 */
enum QUAD_NAME {
    QUAD_MonteCarlo,
    QUAD_GaussLegendreTensorized,
    QUAD_GaussLegendre1D,
    QUAD_GaussLegendreTensorized2D,
    QUAD_LevelSymmetric,
    QUAD_Lebedev,
    QUAD_LDFESA,
    QUAD_Product,
    QUAD_Rectangular1D,
    QUAD_Rectangular2D,
    QUAD_Rectangular3D,
    QUAD_GaussChebyshev1D,
    QUAD_Midpoint1D,
    QUAD_Midpoint2D,
    QUAD_Midpoint3D,
};

/*! @brief Conversion Map String to enum
 */
inline std::map<std::string, QUAD_NAME> Quadrature_Map{ { "MONTE_CARLO", QUAD_MonteCarlo },
                                                        { "GAUSS_LEGENDRE_TENSORIZED", QUAD_GaussLegendreTensorized },
                                                        { "GAUSS_LEGENDRE_TENSORIZED_2D", QUAD_GaussLegendreTensorized2D },
                                                        { "PRODUCT", QUAD_Product },
                                                        { "GAUSS_LEGENDRE_1D", QUAD_GaussLegendre1D },
                                                        { "LEVEL_SYMMETRIC", QUAD_LevelSymmetric },
                                                        { "LEBEDEV", QUAD_Lebedev },
                                                        { "LDFESA", QUAD_LDFESA },
                                                        { "MIDPOINT_1D", QUAD_Midpoint1D },
                                                        { "MIDPOINT_2D", QUAD_Midpoint2D },
                                                        { "MIDPOINT_3D", QUAD_Midpoint3D },
                                                        { "RECTANGULAR_1D", QUAD_Rectangular1D },
                                                        { "RECTANGULAR_2D", QUAD_Rectangular2D },
                                                        { "RECTANGULAR_3D", QUAD_Rectangular3D } };

// Problem name
enum PROBLEM_NAME {
    PROBLEM_LineSource,
    PROBLEM_Checkerboard,
    PROBLEM_ElectronRT,
    PROBLEM_WaterPhantom,
    PROBLEM_LineSource_Pseudo_1D,
    PROBLEM_LineSource_Pseudo_1D_Physics,
    PROBLEM_AirCavity,
    PROBLEM_MuscleBoneLung,
    PROBLEM_Phantom2D,
    PROBLEM_IsotropicSource_2D
};

inline std::map<std::string, PROBLEM_NAME> Problem_Map{ { "LINESOURCE", PROBLEM_LineSource },
                                                        { "CHECKERBOARD", PROBLEM_Checkerboard },
                                                        { "ELECTRONRT", PROBLEM_ElectronRT },
                                                        { "WATERPHANTOM", PROBLEM_WaterPhantom },
                                                        { "AIRCAVITY", PROBLEM_AirCavity },
                                                        { "MUSCLEBONELUNG", PROBLEM_MuscleBoneLung },
                                                        { "PHANTOM2D", PROBLEM_Phantom2D },
                                                        { "ISOTROPICPOINTSOURCE2D", PROBLEM_IsotropicSource_2D },
                                                        { "LINESOURCE_PSEUDO_1D", PROBLEM_LineSource_Pseudo_1D },
                                                        { "LINESOURCE_PSEUDO_1D_PHYSICS", PROBLEM_LineSource_Pseudo_1D_Physics } };

// Kernel name
enum KERNEL_NAME { KERNEL_Isotropic, KERNEL_Isotropic1D };

inline std::map<std::string, KERNEL_NAME> Kernel_Map{ { "ISOTROPIC", KERNEL_Isotropic }, { "ISOTROPIC_1D", KERNEL_Isotropic1D } };

// Solver name
enum SOLVER_NAME {
    SN_SOLVER,
    CSD_SN_SOLVER,
    CSD_SN_NOTRAFO_SOLVER,
    CSD_SN_FOKKERPLANCK_SOLVER,
    CSD_SN_FOKKERPLANCK_TRAFO_SOLVER,
    CSD_SN_FOKKERPLANCK_TRAFO_SOLVER_2D,
    CSD_SN_FOKKERPLANCK_TRAFO_SH_SOLVER_2D,
    CSD_PN_SOLVER,
    PN_SOLVER,
    MN_SOLVER
};

inline std::map<std::string, SOLVER_NAME> Solver_Map{ { "SN_SOLVER", SN_SOLVER },
                                                      { "CSD_SN_SOLVER", CSD_SN_SOLVER },
                                                      { "CSD_SN_NOTRAFO_SOLVER", CSD_SN_NOTRAFO_SOLVER },
                                                      { "CSD_SN_FOKKERPLANCK_SOLVER", CSD_SN_FOKKERPLANCK_SOLVER },
                                                      { "CSD_SN_FOKKERPLANCK_TRAFO_SOLVER", CSD_SN_FOKKERPLANCK_TRAFO_SOLVER },
                                                      { "CSD_SN_FOKKERPLANCK_TRAFO_SOLVER_2D", CSD_SN_FOKKERPLANCK_TRAFO_SOLVER_2D },
                                                      { "CSD_SN_FOKKERPLANCK_TRAFO_SH_SOLVER_2D", CSD_SN_FOKKERPLANCK_TRAFO_SH_SOLVER_2D },
                                                      { "CSD_PN", CSD_PN_SOLVER },
                                                      { "PN_SOLVER", PN_SOLVER },
                                                      { "MN_SOLVER", MN_SOLVER } };

// Entropy functional
enum ENTROPY_NAME { QUADRATIC, MAXWELL_BOLTZMANN, BOSE_EINSTEIN, FERMI_DIRAC };

inline std::map<std::string, ENTROPY_NAME> Entropy_Map{
    { "QUADRATIC", QUADRATIC }, { "MAXWELL_BOLTZMANN", MAXWELL_BOLTZMANN }, { "BOSE_EINSTEIN", BOSE_EINSTEIN }, { "FERMI_DIRAC", FERMI_DIRAC } };

// Optimizer
enum OPTIMIZER_NAME { NEWTON, REGULARIZED_NEWTON, PART_REGULARIZED_NEWTON, ML };

inline std::map<std::string, OPTIMIZER_NAME> Optimizer_Map{
    { "NEWTON", NEWTON }, { "REGULARIZED_NEWTON", REGULARIZED_NEWTON }, { "PART_REGULARIZED_NEWTON", PART_REGULARIZED_NEWTON }, { "ML", ML } };

// Volume output
enum VOLUME_OUTPUT { ANALYTIC, MINIMAL, MOMENTS, DUAL_MOMENTS, MEDICAL };

inline std::map<std::string, VOLUME_OUTPUT> VolOutput_Map{
    { "ANALYTIC", ANALYTIC }, { "MINIMAL", MINIMAL }, { "MOMENTS", MOMENTS }, { "DUAL_MOMENTS", DUAL_MOMENTS }, { "MEDICAL", MEDICAL } };

// Scalar output
enum SCALAR_OUTPUT { ITER, MASS, RMS_FLUX, VTK_OUTPUT, CSV_OUTPUT };

inline std::map<std::string, SCALAR_OUTPUT> ScalarOutput_Map{
    { "ITER", ITER }, { "MASS", MASS }, { "RMS_FLUX", RMS_FLUX }, { "VTK_OUTPUT", VTK_OUTPUT }, { "CSV_OUTPUT", CSV_OUTPUT } };

// Spherical Basis Name
enum SPHERICAL_BASIS_NAME { SPHERICAL_HARMONICS, SPHERICAL_MONOMIALS };

inline std::map<std::string, SPHERICAL_BASIS_NAME> SphericalBasis_Map{ { "SPHERICAL_HARMONICS", SPHERICAL_HARMONICS },
                                                                       { "SPHERICAL_MONOMIALS", SPHERICAL_MONOMIALS } };

// Datasampler Name
enum SAMPLER_NAME { CLASSIFICATION_SAMPLER, REGRESSION_SAMPLER };

inline std::map<std::string, SAMPLER_NAME> SamplerName_MAP{ { "CLASSIFICATION_SAMPLER", CLASSIFICATION_SAMPLER },
                                                            { "REGRESSION_SAMPLER", REGRESSION_SAMPLER } };
#endif    // GLOBAL_CONSTANTS_H
