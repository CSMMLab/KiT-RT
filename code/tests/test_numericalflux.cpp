#include <numeric>

#include "catch.hpp"
#include "fluxes/upwindflux.h"
#include "io.h"
#include "mesh.h"
#include "settings/config.h"
#include "settings/globalconstants.h"
#include "settings/typedef.h"

TEST_CASE( "unit numericalflux tests", "[numericalflux]" ) {
    char config_file_name[] = "../tests/input/unit.cfg";

    Config* config = new Config( config_file_name );

    // Test Flux Jacobians
    Matrix Ax( 4, 4, 0 );
    Matrix AxAbs( 4, 4, 0 );
    Matrix Ay( 4, 4, 0 );
    Matrix AyAbs( 4, 4, 0 );

    // Fill Matrices with used values
    Ax( 0, 3 ) = 0.57735;
    Ax( 3, 0 ) = 0.57735;
    Ay( 0, 1 ) = 0.57735;
    Ay( 1, 0 ) = 0.57735;

    AxAbs( 0, 0 ) = 0.57735;
    AxAbs( 3, 3 ) = 0.57735;
    AyAbs( 0, 0 ) = 0.57735;
    AyAbs( 1, 1 ) = 0.57735;

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

            // ---- Test all fluxes ----

            // a) Upwinding methods - scalar
            UpwindFlux g( config );
            REQUIRE( std::fabs( g.Flux( omega, psiL, psiR, n ) - ( -g.Flux( omega, psiR, psiL, -n ) ) ) <
                     1e2 * std::numeric_limits<double>::epsilon() );

            // b) Upwinding methods - VanLeer  MatrixFlux, 4x4 Matrices, i.e. P1 case;
            Vector resultFluxPlus( 4, 0 );
            Vector resultFluxMinus( 4, 0 );
            Vector psiLVect{ 1.0, 2.0, 3.0, 4.0 };
            Vector psiRVect{ 3.0, 4.0, 1.0, 2.0 };

            g.FluxVanLeer( Ax, AxAbs, Ay, AyAbs, Ay, AyAbs, psiLVect, psiRVect, n, resultFluxPlus );
            g.FluxVanLeer( Ax, AxAbs, Ay, AyAbs, Ay, AyAbs, psiLVect, psiRVect, -n, resultFluxMinus );

            REQUIRE( blaze::norm( resultFluxPlus + resultFluxMinus ) < 1e2 * std::numeric_limits<double>::epsilon() );
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

            // ---- Test all fluxes ----

            // a) Upwinding methods - scalar
            UpwindFlux g( config );

            REQUIRE( std::fabs( g.Flux( omega, psi, psi, n ) - inner * psi ) < 1e2 * std::numeric_limits<double>::epsilon() );

            // b) Upwinding methods - VanLeer  MatrixFlux, 4x4 Matrices, i.e. P1 case;
            Vector resultFlux( 4, 0 );
            Vector resultFluxOrig( 4, 0 );
            Vector psiLVect{ 1.0, 2.0, 3.0, 4.0 };

            resultFluxOrig = ( n[0] * Ax * psiLVect + n[1] * Ay * psiLVect );
            g.FluxVanLeer( Ax, AxAbs, Ay, AyAbs, Ay, AyAbs, psiLVect, psiLVect, n, resultFlux );

            REQUIRE( blaze::norm( resultFlux - resultFluxOrig ) < 1e2 * std::numeric_limits<double>::epsilon() );
        }
    }
}
