#ifndef PROBLEMBASE_H
#define PROBLEMBASE_H

#include "settings/typedef.h"

// Forward Declaration

class Config;
class Physics;
class Mesh;

class ProblemBase
{

  protected:
    Config* _settings;
    Mesh* _mesh;
    Physics* _physics;

    std::vector<double> _density;
    std::vector<double> _stoppingPower;

    ProblemBase() = delete;

  public:
    /**
     * @brief GetScatteringXS gives back vector (each energy) of vectors (each grid cell)
     *        of scattering cross sections for materials defined by density and energies
     *        in vector energy
     * @param energy is the energy the cross section is queried for
     */
    virtual VectorVector GetScatteringXS( const std::vector<double>& energies ) = 0;

    /**
     * @brief GetTotalXS gives back vector of vectors of total cross sections for
     *        materials defined by density and energies in vector energy
     * @param energy is the energy the cross section is queried for
     * @param density is vector with patient densities (at different spatial cells)
     */
    virtual VectorVector GetTotalXS( const std::vector<double>& energies ) = 0;

    /**
     * @brief GetExternalSource gives back vector of vectors of source terms for each
     *        energy, cell and angle
     * @param energies is vector with energies
     */
    virtual std::vector<VectorVector> GetExternalSource( const std::vector<double>& energies ) = 0;

    /**
     * @brief GetStoppingPower gives back vector of vectors of stopping powers for
     *        materials defined by density and energies in vector energy
     * @param energies is vector with energies
     */
    virtual std::vector<double> GetStoppingPower( const std::vector<double>& energies ) = 0;

    /**
     * @brief GetDensity gives back vector of densities for every spatial cell
     * @param cellMidPoints is vector with cell mid points
     */
    virtual std::vector<double> GetDensity( const VectorVector& cellMidPoints );

    /**
     * @brief Setup the initial condition for the flux psi
     */
    virtual VectorVector SetupIC() = 0;

    /**
     * @brief Physics constructor
     * @param settings stores all needed user information
     */
    ProblemBase( Config* settings, Mesh* mesh );
    virtual ~ProblemBase();

    /**
     * @brief Create constructor
     * @param settings stores all needed information
     * @return pointer to Physics
     */
    static ProblemBase* Create( Config* settings, Mesh* mesh );
};

#endif
