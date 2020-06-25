#ifndef PHYSICS_H
#define PHYSICS_H

// include Matrix, Vector definitions
#include "settings/config.h"
#include "settings/typedef.h"
#include <list>
#include <tuple>
#include "math.h"
#include <fstream>
#include <boost/range/combine.hpp>
#include <boost/foreach.hpp>


class Physics
{

  private:
    VectorVector _xsH2O;
    VectorVector _xsTotalH2O;
    VectorVector _xsTransportH2O;
	VectorVector _stpowH2O;

  public:
 
    
	//prototype data readers
    std::tuple<std::list<VectorVector>,std::list<VectorVector>> ReadENDL( std::string filename );
	VectorVector ReadStoppingPowers(std::string fileName);
	
	//load and prepare data from database
	void LoadDatabase (std::string fileName_H, std::string fileName_O, std::string fileName_stppower);


    /** @brief GetScatteringXS gives back vector of vectors of scattering cross sections for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     * @param Omega are scattering angles
     */
    VectorVector GetScatteringXS (Vector energies,Vector density, Vector angle);

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
    VectorVector GetStoppingPower (Vector energies,Vector density);

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
