#ifndef PHYSICS_H
#define PHYSICS_H

// include Matrix, Vector definitions
#include "settings/config.h"
#include "typedef.h"

class Physics
{

  private:
    Matrix _xsH2O;
    Matrix _totalxsH2O;
    /**
     * @brief LoadXSH2O loads (total) scattering cross sections for water from file and saves them in _xsH2O and _totalxsH2O
     * @param fileName1 is name of cross section file
     * @param fileName2 is name of total cross section file
     */
    void LoadXSH2O( std::string fileName1, std::string fileName2 );

  public:
    /**
     * @brief GetScatteringXS gives back vector of vectors of scattering cross sections for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     * @param Omega are scattering angles
     */
    VectorVector GetScatteringXS( std::vector<double> energies, std::vector<double> density, std::vector<double> Omegas );

    /**
     * @brief GetTotalXS gives back vector of vectors of total cross sections for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     */
    VectorVector GetTotalXS( std::vector<double> energies, std::vector<double> density );

    /**
     * @brief GetStoppingPower gives back vector of vectors of stopping powers for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     * @param sH2O is vector of stopping powers in water
     */
    VectorVector GetStoppingPower( std::vector<double> energies, std::vector<double> sH2O );

    /**
     * @brief Physics constructor
     * @param settings stores all needed user information
     */
    Physics( Config* settings );

    /**
     * @brief Create constructor
     * @param settings stores all needed information
     * @return pointer to Physics
     */
    static Physics* Create( Config* settings );
};

#endif
