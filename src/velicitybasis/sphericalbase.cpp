/*!
 * @file sphericalbase.cpp
 * @brief Base Class to handle basis classes on the unit sphere
 * @author S. Schotthöfer
 */

#include "velocitybasis/sphericalbase.hpp"
#include "common/config.hpp"
#include "toolboxes/errormessages.hpp"
#include "velocitybasis/sphericalharmonics.hpp"
#include "velocitybasis/sphericalmonomials.hpp"

SphericalBase* SphericalBase::Create( Config* settings ) {
    SPHERICAL_BASIS_NAME name = settings->GetSphericalBasisName();
    unsigned maxMomentDegree  = settings->GetMaxMomentDegree();
    unsigned short spatialDim = settings->GetDim();

    switch( name ) {
        case SPHERICAL_HARMONICS: return new SphericalHarmonics( maxMomentDegree, spatialDim ); break;
        case SPHERICAL_MONOMIALS: return new SphericalMonomials( maxMomentDegree, spatialDim ); break;
        default: ErrorMessages::Error( "Creator for the chosen basis does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    }
    return nullptr;
}
