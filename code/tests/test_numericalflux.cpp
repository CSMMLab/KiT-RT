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
    Matrix AxP( 4, 4, 0 );
    Matrix AxM( 4, 4, 0 );
    Matrix AyP( 4, 4, 0 );
    Matrix AyM( 4, 4, 0 );

    // Fill Matrices with used values
    AxP( 0, 0 ) = 1;    // 0.5;
    AxP( 1, 1 ) = 0;    // 0.25;
    AyP( 0, 0 ) = 0;    // 0.25;
    AyP( 1, 1 ) = 1;    // 0.125;

    AxM( 2, 2 ) = -1;    // -0.25;
    AxM( 3, 3 ) = 0;     // -0.125;
    AyM( 2, 2 ) = 0;     // -0.5;
    AyM( 3, 3 ) = -1;    // -0.25;

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
        }
        // b) Upwinding methods - Systems  MatrixFlux, 4x4 Matrices, i.e. P1 case;
        Vector resultFluxPlus( 4, 0.0 );
        Vector resultFluxMinus( 4, 0.0 );
        Vector psiLVect{ 1.0, 2.0, 3.0, 4.0 };
        Vector psiRVect{ 3.0, 4.0, 1.0, 2.0 };
        // n[0] = 1.0;
        // n[1] = 1.0;
        // std::cout << "resultFluxPlus " << resultFluxPlus << " \n resultFluxMinus " << resultFluxMinus << "\n";
        // std::cout << " n " << n << "\n";

        // std::cout << "AxP \n" << AxP << "\n AxM \n" << AxM << "\n AyP \n" << AyP << "\n AyM \n" << AyM << "\n";
        resultFluxPlus  = g.Flux( AxP, AxM, AyP, AyM, AyP, AyM, psiLVect, psiRVect, n );
        resultFluxMinus = g.Flux( AxP, AxM, AyP, AyM, AyP, AyM, psiRVect, psiLVect, -n );
        // std::cout << "resultFluxPlus " << resultFluxPlus << " \n resultFluxMinus " << resultFluxMinus << "\n";

        REQUIRE( blaze::norm( resultFluxPlus + resultFluxMinus ) < 1e2 * std::numeric_limits<double>::epsilon() );
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
        }
        // b) Upwinding methods -  MatrixFlux, 4x4 Matrices, i.e. P1 case;
        Vector resultFlux( 4, 0 );
        Vector resultFluxOrig( 4, 0 );
        Vector psiLVect{ 1.0, 2.0, 3.0, 4.0 };

        // std::cout << "resultFluxPlus " << resultFlux << "resOrig " << resultFluxOrig << "\n";
        // std::cout << " n " << n << "\n";
        // std::cout << "Ax \n" << AxP + AxM << "\n Ay \n" << AyP + AyM << "\n";

        resultFluxOrig = ( n[0] * ( AxP + AxM ) * psiLVect + +n[1] * ( AyP + AyM ) * psiLVect );
        resultFlux     = g.Flux( AxP, AxM, AyP, AyM, AyP, AyM, psiLVect, psiLVect, n );

        REQUIRE( blaze::norm( resultFlux - resultFluxOrig ) < 1e2 * std::numeric_limits<double>::epsilon() );
    }
}
