#ifndef PHYSICS_H
#define PHYSICS_H

#include <fstream>
#include <list>
#include <map>

#include "common/typedef.h"
#include "spline.h"
#include "toolboxes/errormessages.h"

class Physics
{

  private:
    VectorVector _xsH2O;
    VectorVector _xsTotalH2O;
    VectorVector _xsTransportH2O;
    VectorVector _stpowH2O;

    std::tuple<std::vector<VectorVector>, std::vector<VectorVector>> ReadENDL( std::string filename );
    VectorVector ReadStoppingPowers( std::string fileName );
    void LoadDatabase( std::string fileName_H, std::string fileName_O, std::string fileName_stppower );

    Physics() = delete;

  public:
    /** @brief GetScatteringXS gives back vector of vectors of scattering cross sections for materials defined by density and energies in vector
     * energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     * @param Omega are scattering angles
     */
    VectorVector GetScatteringXS( Vector energies, Vector angle );

    VectorVector GetScatteringXSE( Vector energies );

    /**
     * @brief GetTotalXS gives back vector of vectors of total cross sections for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     */
    VectorVector GetTotalXS( Vector energies, Vector density );

    VectorVector GetTotalXSE( Vector energies );

    /**
     * @brief GetStoppingPower gives back vector of vectors of stopping powers for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     * @param sH2O is vector of stopping powers in water
     */
    Vector GetStoppingPower( Vector energies );

    /**
     * @brief GetTransportXS gives back vector of vectors of stopping powers for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     */
    VectorVector GetTransportXS( Vector energies, Vector density );

    VectorVector GetTransportXSE( Vector energies );

    /**
     * @brief Physics constructor
     */
    Physics( std::string fileName_H, std::string fileName_O, std::string fileName_stppower );

    /**
     * @brief Physics destructor
     */
    ~Physics();
};

#endif
