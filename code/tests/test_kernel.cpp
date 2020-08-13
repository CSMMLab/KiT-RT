#include <numeric>

#include "catch.hpp"
#include "common/config.h"
#include "kernels/isotropic.h"
#include "quadratures/quadraturebase.h"

TEST_CASE( "test all scattering kernels", "[kernel]" ) {
    char filename[] = "../tests/input/unit_kernel.cfg";

    // Load Settings from File
    Config* config = new Config( filename );

    QuadratureBase* quad = QuadratureBase::CreateQuadrature( config );    //@TODO: swap out for different quadrature rule

    SECTION( "isotropic scattering kernel" ) {

        auto weights = quad->GetWeights();
        Isotropic kernel( quad );
        Matrix scatteringMatrix = kernel.GetScatteringKernel();
        for( unsigned i = 0; i < scatteringMatrix.rows(); ++i ) {
            for( unsigned j = 0; j < scatteringMatrix.columns(); ++j ) {
                REQUIRE( std::fabs( scatteringMatrix( i, j ) - ( weights[j] / ( 4 * M_PI ) ) ) < std::numeric_limits<double>::epsilon() );
            }
        }
    }
}
