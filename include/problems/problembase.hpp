#ifndef PROBLEMBASE_H
#define PROBLEMBASE_H

#include "common/typedef.hpp"

// Forward Declaration

class Config;
class EPICS;
class Mesh;

class ProblemBase
{

  protected:
    Config* _settings;
    Mesh* _mesh;
    EPICS* _physics;

    std::vector<double> _density;       /*!< @brief vector with patient densities */
    std::vector<double> _stoppingPower; /*!< @brief vector with stopping powers*/

    ProblemBase() = delete;

  public:
    /**
     * @brief GetScatteringXS gives back vector (each energy) of vectors (each grid cell)
     *        of scattering cross sections for materials defined by density and energies
     *        in vector energy
     * @param energies is the energy the cross section is queried for
     */
    virtual VectorVector GetScatteringXS( const Vector& energies ) = 0;

    /**
     * @brief GetTotalXS gives back vector of vectors of total cross sections for
     *        materials defined by density and energies in vector energy
     * @param energies is the energy the cross section is queried for
     */
    virtual VectorVector GetTotalXS( const Vector& energies ) = 0;

    /**
     * @brief GetTotalXSE gives back vector of total cross sections for
     *        energies in vector energy
     */
    virtual Vector GetTotalXSE( const Vector& /*energies*/ ) { return Vector( 1 ); }

    /**
     * @brief GetScatteringXSE gives back vector (each energy) of scattering cross sections for energies
     *        in vector energy
     * @param energies is the energy the cross section is queried for
     * @param angles are the queried angles
     */
    virtual std::vector<Matrix> GetScatteringXSE( const Vector& energies, const Matrix& angles ) {
        return std::vector<Matrix>( energies.size(), Matrix( angles.rows(), angles.columns() ) );
    }

    /**
     * @brief GetScatteringXSE gives back vector (each energy) of scattering cross sections for energies
     *        in vector energy
     * @param energies is the energy the cross section is queried for
     * @param angles are the queried angles
     */
    virtual VectorVector GetScatteringXSE( const Vector& energies, const Vector& angles );

    /**
     * @brief GetExternalSource gives back vector of vectors of source terms for each
     *        energy, cell and angle
     * @param energies is the energy the cross section is queried for
     */
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) = 0;

    /**
     * @brief GetStoppingPower gives back vector of vectors of stopping powers for
     *        materials defined by density and energies in vector energy
     * @param energies is vector with energies
     */
    virtual Vector GetStoppingPower( const Vector& energies );

    /**
     * @brief GetDensity gives back vector of densities for every spatial cell
     * @param cellMidPoints is vector with cell mid points
     */
    virtual std::vector<double> GetDensity( const VectorVector& cellMidPoints );

    /**
     * @brief Setup the initial condition for the flux psi
     */
    virtual VectorVector SetupIC() = 0;

    /*! @brief Exact analytical solution for the Line Source Test Case. Returns 0 for all other test cases.
         @return exact solution at x,y,t,scatteringXS
    */ // Default is set to 0. ~> if no analytical solution is available.
    double virtual GetAnalyticalSolution( double /*x*/, double /*y*/, double /*t*/, double /*scatteringXS*/ ) { return 0.0; }

    /**
     * @brief Physics constructor
     * @param settings stores all needed user information
     * @param mesh for the test case
     */
    ProblemBase( Config* settings, Mesh* mesh );
    virtual ~ProblemBase();

    /**
     * @brief Create constructor
     * @param settings stores all needed information
     * @param mesh for the test case
     * @return pointer to ProblemBase
     */
    static ProblemBase* Create( Config* settings, Mesh* mesh );
};

#endif
