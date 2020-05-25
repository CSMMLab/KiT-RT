#ifndef SOLVER_H
#define SOLVER_H

#include <string>

// include Matrix, Vector definitions
#include "io.h"
#include "kernels/scatteringkernelbase.h"
#include "numericalflux.h"
#include "problems/problembase.h"
#include "quadratures/quadraturebase.h"
#include "settings/config.h"
#include "typedef.h"

class Solver
{
  protected:
    unsigned _nq;                                     // number of quadrature points
    unsigned _nCells;                                 // number of spatial cells
    unsigned _nEnergies;                              // number of energy/time steps, number of nodal energy values for CSD
    double _dE;                                       // energy/time step size
    std::vector<double> _energies;                    // energy groups used in the simulation [keV]
    VectorVector _psi;                                // angular flux vector, dim(_psi) = (_NCells,_nq)
    std::vector<double> _areas;                       // surface area of all spatial cells, dim(_areas) = _NCells
    std::vector<std::vector<Vector>> _normals;        // edge normals multiplied by edge length, dim(_normals) = (_NCells,nEdgesPerCell,spatialDim)
    std::vector<std::vector<unsigned>> _neighbors;    // edge normals multiplied by edge length, dim(_neighbors) = (_NCells,nEdgesPerCell)
    std::vector<double> _density;                     // patient density, dim(_density) = _nCells
    std::vector<double> _s;                           // stopping power, dim(_s) = _nTimeSteps
    VectorVector _sigmaS;                             // scattering cross section for all energies
    VectorVector _sigmaT;                             // total cross section for all energies
    VectorVector _Q;                                  // external source term
    Matrix _scatteringKernel;                         // scattering kernel for the quadrature
    VectorVector _quadPoints;                         // quadrature points, dim(_quadPoints) = (_nTimeSteps,spatialDim)
    Vector _weights;                                  // quadrature weights, dim(_weights) = (_NCells)
    std::vector<BOUNDARY_TYPE> _boundaryCells;        // boundary type for all cells, dim(_boundary) = (_NCells)
    // we will have to add a further dimension for quadPoints and weights once we start with multilevel SN

    NumericalFlux* _g;    // numerical flux function
    Mesh* _mesh;          // mesh object for writing out information
    Config* _settings;
    ProblemBase* _problem;

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
};

#endif    // SOLVER_H
