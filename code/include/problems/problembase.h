#ifndef PROBLEMBASE_H
#define PROBLEMBASE_H

#include "mesh.h"
#include "physics.h"
#include "settings/config.h"
#include "typedef.h"

class ProblemBase
{

  protected:
    Config* _settings;
    Mesh* _mesh;
    Physics* _physics;

    std::vector<double> _totalXS;
    Matrix _scatteringXS;
    std::vector<double> _density;
    std::vector<double> _stoppingPower;

    ProblemBase() = delete;

  public:
    /**
     * @brief GetScatteringXS gives back vector of vectors of scattering cross sections for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     * @param Omega are scattering angles
     */
    virtual std::vector<Matrix> GetScatteringXS( const double energy ) = 0;

    /**
     * @brief GetTotalXS gives back vector of vectors of total cross sections for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     */
    virtual std::vector<double> GetTotalXS( const double energy ) = 0;

    /**
     * @brief GetStoppingPower gives back vector of vectors of stopping powers for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     * @param sH2O is vector of stopping powers in water
     */
    virtual std::vector<double> GetStoppingPower( const std::vector<double>& energies ) = 0;

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
