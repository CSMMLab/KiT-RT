/*!
 * @file sphericalbase.cpp
 * @brief Base Class to handle basis classes on the unit sphere
 * @author S. Schotthöfer
 */

#include "toolboxes/sphericalbase.h"
#include "common/config.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/sphericalharmonics.h"
#include "toolboxes/sphericalmonomials.h"

SphericalBase* SphericalBase::Create( Config* settings ) {
    SPHERICAL_BASIS_NAME name = settings->GetSphericalBasisName();
    unsigned maxMomentDegree  = settings->GetMaxMomentDegree();
    unsigned short spatialDim = settings->GetDim();

    switch( name ) {
        case SPHERICAL_HARMONICS:
            if( spatialDim == 3 ) {
                return new SphericalHarmonics( maxMomentDegree );
            }
            ErrorMessages::Error( "Spherical Harmonics basis is not yet equipped for 1D and 2D cases.", CURRENT_FUNCTION );
            break;
        case SPHERICAL_MONOMIALS: return new SphericalMonomials( maxMomentDegree, spatialDim );

        default: ErrorMessages::Error( "Creator for the chosen basis does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    }
    return nullptr;
}
