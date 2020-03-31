#ifndef SOLVER_H
#define SOLVER_H

// include Matrix, Vector definitions
#include "typedef.h"

#include "numericalflux.h"
#include "quadrature.h"
#include "settings.h"
#include <string>

class Solver
{
  protected:
    unsigned _nq;                                     // number of quadrature points
    const unsigned _NCells;                           // number of spatial cells
    const unsigned _NTimeSteps;                       // number of time steps, number of nodal energy values for CSD
    double _dt;                                       // time step size
    Matrix _psi;                                      // angular flux vector, dim(_psi) = (_NCells,_nq)
    std::vector<unsigned> _areas;                     // surface area of all spatial cells, dim(_areas) = _NCells
    std::vector<std::vector<Vector>> _normals;        // edge normals multiplied by edge length, dim(_normals) = (_NCells,spatialDim)
    std::vector<std::vector<unsigned>> _neighbors;    // edge normals multiplied by edge length, dim(_neighbors) = (_NCells,nEdgesPerCell)
    std::vector<double> _sH20;                        // stopping power H2O, dim(_sH20) = _nTimeSteps
    std::vector<double> _density;                     // patient density, dim(_density) = _nCells
    VectorVector _quadPoints;                         // quadrature points, dim(_quadPoints) = (_NCells,spatialDim)
    Vector _weights;                                  // quadrature weights, dim(_weights) = (_NCells)
    // we will have to add a further dimension for quadPoints and weights once we start with multilevel SN

    NumericalFlux* _g;    // numerical flux function

    /**
     * @brief LoadPatientDensity loads density of patient from MRT/CT scan and saves it in _density
     * @param fileName is name of patient file
     */
    void LoadPatientDensity( std::string fileName );

    /**
     * @brief LoadStoppingPower loads stopping power of H20 and saves it in _sH20
     * @param fileName is name of stopping power file
     */
    void LoadStoppingPower( std::string fileName );

  public:
    /**
     * @brief Solver constructor
     * @param settings stores all needed information
     */
    Solver( Settings* settings );

    /**
     * @brief Create constructor
     * @param settings stores all needed information
     * @return pointer to Solver
     */
    static Solver* Create( Settings* settings );

    /**
     * @brief Solve functions runs main time loop
     */
    virtual void Solve() = 0;
};

#endif    // SOLVER_H
