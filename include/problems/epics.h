#ifndef EPICS_H
#define EPICS_H

#include <fstream>
#include <list>
#include <map>

#include "common/globalconstants.h"
#include "common/typedef.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/interpolation.h"

class EPICS
{
  private:
    enum Element { H = 0, O = 1 };
    std::vector<VectorVector> _xsScatteringH2O;    // only holds angular distribution, no cross section terms
    std::vector<VectorVector> _xsTotalH2O;         // equal to scatteringXS integrated over angle
    VectorVector _stpowH2O;

    std::tuple<std::vector<VectorVector>, std::vector<VectorVector>> ReadEPICS( std::string filename );
    VectorVector ReadStoppingPowers( std::string fileName );
    void LoadDatabase( std::string fileName_H, std::string fileName_O, std::string fileName_stppower );
    VectorVector ReorderTotalXS( const VectorVector& data ) const;
    VectorVector ReorderScatteringXS( const VectorVector& data ) const;

    EPICS() = delete;

  public:
    /** @brief GetScatteringXS gives back vector of vectors of scattering cross sections for materials defined by density and energies in vector
     * energy
     * @param energies is vector with energies
     * @param angle are scattering angles (mu)
     */
    VectorVector GetScatteringXS( Vector energies, Vector angle );

    /** @brief GetScatteringXS gives back vector of Matrix of scattering cross sections for materials defined by density and energies in vector
     * energy
     * @param energies is vector with energies
     * @param angle is matrix of scattering angles (<Omega_i,Omega_j>)
     */
    std::vector<Matrix> GetScatteringXS( const Vector& energies, const Matrix& angle );

    /**
     * @brief GetTotalXS gives back vector of vectors of total cross sections for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     * @param density is vector with patient densities (at different spatial cells)
     */
    VectorVector GetTotalXS( Vector energies, Vector density );

    Vector GetTotalXSE( Vector energies );

    /**
     * @brief GetStoppingPower gives back vector of vectors of stopping powers for materials defined by density and energies in vector energy
     * @param energies is vector with energies
     */
    Vector GetStoppingPower( Vector energies );

    /**
     * @brief EPICS constructor
     */
    EPICS( std::string fileName_H, std::string fileName_O, std::string fileName_stppower = "" );

    /**
     * @brief EPICS destructor
     */
    ~EPICS();
};

#endif
