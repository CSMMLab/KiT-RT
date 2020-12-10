/*!
 * @file sphericalbasisbase.cpp
 * @brief Base Class to handle basis classes on the unit sphere
 * @author S. Schotthöfer
 *
 */

#include "toolboxes/sphericalbasisbase.h"
#include "common/config.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/sphericalharmonics.h"
#include "toolboxes/sphericalmonomials.h"

SphericalBasisBase* SphericalBasisBase::Create( Config* settings ) {
    SPHERICAL_BASIS_NAME name = settings->GetSphericalBasisName();
    unsigned maxMomentDegree  = settings->GetMaxMomentDegree();

    switch( name ) {
        case SPHERICAL_HARMONICS: return new SphericalHarmonics( maxMomentDegree );
        case SPHERICAL_MONIMIALS: return new SphericalMonomials( maxMomentDegree );

        default: ErrorMessages::Error( "Creator for the chosen basis does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    }
    return nullptr;
}
