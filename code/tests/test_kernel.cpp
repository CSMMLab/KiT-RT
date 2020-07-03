#include <numeric>

#include "catch.hpp"
#include "kernels/isotropic.h"
#include "quadratures/qlebedev.h"

TEST_CASE( "test all scattering kernels", "[kernel]" ) {
    unsigned nq    = 5;
    QLebedev* quad = new QLebedev( nq );    //@TODO: swap out for different quadrature rule

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
