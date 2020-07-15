#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
// include Matrix, Vector definitions
#include "settings/globalconstants.h"
#include "settings/typedef.h"

// Forward Declarations
class NumericalFlux;
class Mesh;
class Config;
class ProblemBase;
class QuadratureBase;

class Solver
{
  protected:
    Mesh* _mesh;           /*! @brief mesh object for writing out information */
    NumericalFlux* _g;     /*! @brief class for numerical flux */
    Config* _settings;     /*! @brief config class for global information */
    ProblemBase* _problem; /*! @brief problem class for initial conditions */

    // --------- Often used variables of member classes for faster access ----

    unsigned _nEnergies;              /*! @brief number of energy/time steps, number of nodal energy values for CSD */
    double _dE;                       /*! @brief energy/time step size */
    std::vector<double> _energies;    // energy groups used in the simulation [keV]
    std::vector<double> _density;     // patient density, dim(_density) = _nCells
    std::vector<double> _s;           // stopping power, dim(_s) = _nTimeSteps
    std::vector<VectorVector> _Q;     /*!  @brief  external source term */

    VectorVector _sigmaS; /*!  @brief scattering cross section for all energies */
    VectorVector _sigmaT; /*!  @brief total cross section for all energies */

    // quadrature related numbers
    QuadratureBase* _quadrature; /*! @brief quadrature to create members below */
    unsigned _nq;                /*! @brief number of quadrature points */

    // VectorVector _quadPoints;    /*!  @brief quadrature points, dim(_quadPoints) = (_nSystem,spatialDim) */
    // Vector _weights;             /*!  @brief quadrature weights, dim(_weights) = (_NCells) */

    // Mesh related members
    unsigned _nCells;                          /*! @brief number of spatial cells */
    std::vector<BOUNDARY_TYPE> _boundaryCells; /*! boundary type for all cells, dim(_boundary) = (_NCells) */
    std::vector<double> _areas;                /*! @brief surface area of all spatial cells, dim(_areas) = _NCells */
    /*! @brief edge normals multiplied by edge length, dim(_normals) = (_NCells,nEdgesPerCell,spatialDim) */
    std::vector<std::vector<Vector>> _normals;
    /*! @brief edge neighbor cell ids, dim(_neighbors) = (_NCells,nEdgesPerCell) */
    std::vector<std::vector<unsigned>> _neighbors;

    // Solution related members
    VectorVector _sol;                 /*! @brief solution of the PDE, e.g. angular flux or moments */
    std::vector<double> _solverOutput; /*! @brief PROTOTYPE: Outputfield for solver */

    // we will have to add a further dimension for quadPoints and weights once we start with multilevel SN

    /**
     * @brief ComputeTimeStep calculates the maximal stable time step
     * @param cfl is cfl number
     */
    double ComputeTimeStep( double cfl ) const;

  public:
    /**
     * @brief Solver constructor
     * @param settings stores all needed information
     */
    Solver( Config* settings );

    ~Solver();

    /**
     * @brief Create constructor
     * @param settings stores all needed information
     * @return pointer to Solver
     */
    static Solver* Create( Config* settings );

    /**
     * @brief Solve functions runs main time loop
     */
    virtual void Solve() = 0;

    /**
     * @brief Output solution to VTK file
     */
    virtual void Save() const = 0;

    virtual void Save( int currEnergy ) const = 0;
};

#endif    // SOLVER_H
