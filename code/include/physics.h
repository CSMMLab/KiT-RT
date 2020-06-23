#ifndef PHYSICS_H
#define PHYSICS_H

// include Matrix, Vector definitions
#include "settings/config.h"
#include "settings/typedef.h"

#include "math.h"
#include <fstream>


class Physics
{

  private:
    VectorVector _xsH2O;
    VectorVector _xsTotalH2O;
    VectorVector _xsTransportH2O;
	VectorVector _stpowH2O;
    /**
     * @brief LoadDatabase loads required physics data from file in ENDL format (default IAEA EEDL database)
     * @param fileName is name of cross section file
     */

  public:
 
    void LoadDatabase (std::string fileName);
    
     /**
     * @brief GetScatteringXS gives back vector of vectors of scattering cross sections for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     * @param Omega are scattering angles
     */
    VectorVector GetScatteringXS (Vector energies,Vector density, Vector Omegas);

    /**
     * @brief GetTotalXS gives back vector of vectors of total cross sections for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     */
    VectorVector GetTotalXS (Vector energies,Vector density);

    /**
     * @brief GetStoppingPower gives back vector of vectors of stopping powers for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     * @param sH2O is vector of stopping powers in water
     */
    VectorVector GetStoppingPower (Vector energies,Vector density,Vector sH2O);

    /**
     * @brief GetTransportXS gives back vector of vectors of stopping powers for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     */
    VectorVector GetTransportXS (Vector energies,Vector density);

    /**
     * @brief Physics constructor
     * @param settings stores all needed user information
     */
    Physics ();

    /**
     * @brief Create constructor
     * @param settings stores all needed information
     * @return pointer to Physics
     */
    static Physics* Create ();
};

#endif
