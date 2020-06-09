#include <numeric>

#include "catch.hpp"
#include "io.h"
#include "mesh.h"
#include "numericalflux.h"
#include "settings/config.h"
#include "settings/globalconstants.h"

TEST_CASE( "unit numericalflux tests", "[numericalflux]" ) {
    char config_file_name[] = "../tests/input/unit.cfg";

    Config* config = new Config( config_file_name );

    // setup numerical flux
    NumericalFlux* g = NumericalFlux::Create( config );

    SECTION( "test symmetry" ) {
        for( unsigned l = 0; l < 10; ++l ) {
            Vector n( 2, 0.0 );
            std::default_random_engine generator;
            std::normal_distribution<double> distribution( -10.0, 10.0 );
            // sample normal
            n[0] = distribution( generator );
            n[1] = distribution( generator );
            // sample edge length
            double length = distribution( generator );

            // scale normalized normal
            n = length * n / sqrt( n[0] * n[0] + n[1] * n[1] );

            // sample solution at neighboring cells
            double psiL = distribution( generator ) + 10.0;
            double psiR = distribution( generator ) + 10.0;

            // sample omega (note that omega now does not need to lie on unit sphere)
            Vector omega( 2, 0.0 );
            std::normal_distribution<double> uni( -1.0, 1.0 );
            omega[0] = uni( generator );
            omega[1] = uni( generator );

            REQUIRE( std::fabs( g->Flux( omega, psiL, psiR, n ) - ( -g->Flux( omega, psiR, psiL, -n ) ) ) <
                     1e2 * std::numeric_limits<double>::epsilon() );
        }
    }

    SECTION( "test consistency" ) {
        for( unsigned l = 0; l < 10; ++l ) {
            Vector n( 2, 0.0 );
            std::default_random_engine generator;
            std::normal_distribution<double> distribution( -10.0, 10.0 );
            // sample normal
            n[0] = distribution( generator );
            n[1] = distribution( generator );
            // sample edge length
            double length = distribution( generator );

            // scale normalized normal
            n = length * n / sqrt( n[0] * n[0] + n[1] * n[1] );

            // sample solution at neighboring cells
            double psi = distribution( generator ) + 10.0;

            // sample omega (note that omega now does not need to lie on unit sphere)
            Vector omega( 2, 0.0 );
            std::normal_distribution<double> uni( -1.0, 1.0 );
            omega[0] = uni( generator );
            omega[1] = uni( generator );

            double inner = omega[0] * n[0] + omega[1] * n[1];

            REQUIRE( std::fabs( g->Flux( omega, psi, psi, n ) - inner * psi ) < 1e2 * std::numeric_limits<double>::epsilon() );
        }
    }
}
