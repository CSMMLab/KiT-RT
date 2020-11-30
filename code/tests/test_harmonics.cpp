
#include "catch.hpp"
#include "common/config.h"
#include "quadratures/qgausslegendretensorized.h"
#include "solvers/sphericalharmonics.h"

#include <fstream>
#include <iostream>
#include <sstream>

double Y0_0( double, double ) { return sqrt( 1 / ( 4 * M_PI ) ); }

double Y1_m1( double my, double phi ) { return -sqrt( 3 / ( 4 * M_PI ) ) * sqrt( 1 - my * my ) * sin( phi ); }
double Y1_0( double my, double /*phi*/ ) { return sqrt( 3 / ( 4 * M_PI ) ) * my; }
double Y1_1( double my, double phi ) { return -sqrt( 3 / ( 4 * M_PI ) ) * sqrt( 1 - my * my ) * cos( phi ); }

double Y2_m2( double my, double phi ) { return sqrt( 15 / ( 16 * M_PI ) ) * ( 1 - my * my ) * sin( 2 * phi ); }
double Y2_m1( double my, double phi ) { return -1 * sqrt( 15 / ( 4 * M_PI ) ) * my * sqrt( 1 - my * my ) * sin( phi ); }
double Y2_0( double my, double /*phi*/ ) { return sqrt( 5 / ( 16 * M_PI ) ) * ( 3 * my * my - 1 ); }
double Y2_1( double my, double phi ) { return -1 * sqrt( 15 / ( 4 * M_PI ) ) * my * sqrt( 1 - my * my ) * cos( phi ); }
double Y2_2( double my, double phi ) { return sqrt( 15 / ( 16 * M_PI ) ) * ( 1 - my * my ) * cos( 2 * phi ); }

double P0_0( double /*my*/ ) { return sqrt( 1 / ( 2 * M_PI ) ); }
double P1_0( double my ) { return sqrt( 3 / ( 2 * M_PI ) ) * my; }
double P1_1( double my ) { return -sqrt( 3 / ( 4 * M_PI ) ) * sqrt( 1 - my * my ); }

double P2_0( double my ) { return sqrt( 5 / ( 8 * M_PI ) ) * ( 3 * my * my - 1 ); }
double P2_1( double my ) { return -1 * sqrt( 15 / ( 4 * M_PI ) ) * my * sqrt( 1 - my * my ); }
double P2_2( double my ) { return sqrt( 15 / ( 16 * M_PI ) ) * ( 1 - my * my ); }

TEST_CASE( "test the spherical harmonics basis computation", "[spherical_harmonics]" ) {

    std::string filename = std::string( TESTS_PATH ) + "input/unit_tests/solvers/unit_harmonics.cfg";

    // Load Settings from File
    Config* config = new Config( filename );

    unsigned maxMomentDegree = 2;

    SphericalHarmonics testBase( maxMomentDegree );

    SECTION( "Test against analytical solution" ) {
        std::vector<double> legendre;
        Vector moment;
        std::vector<bool> validLegendrePoly( 6, true );
        std::vector<bool> validMoment( 9, true );
        for( double my = -1.0; my < 1.0; my += 0.1 ) {

            legendre = testBase.GetAssLegendrePoly( my );

            if( std::fabs( legendre[0] - P0_0( my ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validLegendrePoly[0] = false;
            if( std::fabs( legendre[1] - P1_0( my ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validLegendrePoly[1] = false;
            if( std::fabs( legendre[2] - P1_1( my ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validLegendrePoly[2] = false;
            if( std::fabs( legendre[3] - P2_0( my ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validLegendrePoly[3] = false;
            if( std::fabs( legendre[4] - P2_1( my ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validLegendrePoly[4] = false;
            if( std::fabs( legendre[5] - P2_2( my ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validLegendrePoly[5] = false;

            for( double phi = 0.0; phi < 2 * M_PI; phi += 0.1 ) {
                moment = testBase.ComputeSphericalBasis( my, phi );

                if( std::fabs( moment[0] - Y0_0( my, phi ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validMoment[0] = false;
                if( std::fabs( moment[1] - Y1_m1( my, phi ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validMoment[1] = false;
                if( std::fabs( moment[2] - Y1_0( my, phi ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validMoment[2] = false;
                if( std::fabs( moment[3] - Y1_1( my, phi ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validMoment[3] = false;
                if( std::fabs( moment[4] - Y2_m2( my, phi ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validMoment[4] = false;
                if( std::fabs( moment[5] - Y2_m1( my, phi ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validMoment[5] = false;
                if( std::fabs( moment[6] - Y2_0( my, phi ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validMoment[6] = false;
                if( std::fabs( moment[7] - Y2_1( my, phi ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validMoment[7] = false;
                if( std::fabs( moment[8] - Y2_2( my, phi ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) validMoment[8] = false;
            }
        }
        REQUIRE( std::all_of( validLegendrePoly.begin(), validLegendrePoly.end(), []( bool v ) { return v; } ) );
        REQUIRE( std::all_of( validMoment.begin(), validMoment.end(), []( bool v ) { return v; } ) );
    }

    // Remove title line

    SECTION( "test to reference solution" ) {

        std::string text_line;
        std::ifstream case_file;

        double my  = 0.0;
        double phi = 0.0;

        Vector values( 9, 0.0 );
        Vector result( 4, 0.0 );

        case_file.open( "unit_test/solvers/harmonicBasis_reference.csv", std::ios::in );

        getline( case_file, text_line );

        bool errorWithinBounds = true;
        while( getline( case_file, text_line ) ) {

            // give line to stringstream
            std::stringstream ss( text_line );

            // Read values
            ss >> my >> phi >> values[0] >> values[1] >> values[2] >> values[3] >> values[4] >> values[5] >> values[6] >> values[7] >> values[8];

            result = testBase.ComputeSphericalBasis( my, phi );

            for( unsigned idx = 0; idx < 9; idx++ ) {
                if( std::fabs( result[idx] - values[idx] ) > 1e2 * std::numeric_limits<double>::epsilon() ) errorWithinBounds = false;
            }
        }
        REQUIRE( errorWithinBounds );
        case_file.close();
    }

    SECTION( "test orthonormality - spherical coordinates" ) {
        // Caution: Integration only works with spherical coordinates!

        QGaussLegendreTensorized quad( config );

        double my, phi, w;
        Vector moment = testBase.ComputeSphericalBasis( 0, 1, 0 );
        // 9 basis moments if degree = 2

        Matrix results( moment.size(), moment.size(), 0.0 );

        for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
            my  = quad.GetPointsSphere()[idx_quad][0];
            phi = quad.GetPointsSphere()[idx_quad][1];
            // z      = quad.GetPoints()[idx_quad][2];
            w      = quad.GetWeights()[idx_quad];
            moment = testBase.ComputeSphericalBasis( my, phi );

            for( unsigned idx_row = 0; idx_row < 9; idx_row++ ) {
                for( unsigned idx_col = 0; idx_col < 9; idx_col++ ) {
                    results( idx_row, idx_col ) += w * moment[idx_row] * moment[idx_col];
                }
            }
        }

        bool errorWithinBounds = true;
        bool orthogonality     = true;
        for( unsigned idx_row = 0; idx_row < 9; idx_row++ ) {
            for( unsigned idx_col = 0; idx_col < 9; idx_col++ ) {
                if( idx_row == idx_col ) {
                    // Orthogonality
                    if( std::fabs( results( idx_row, idx_col ) - 1.0 ) > 1e2 * std::numeric_limits<double>::epsilon() ) orthogonality = false;
                }
                else {
                    // Normality
                    if( std::fabs( results( idx_row, idx_col ) ) > 1e2 * std::numeric_limits<double>::epsilon() ) errorWithinBounds = false;
                }
            }
        }
        REQUIRE( errorWithinBounds );
        REQUIRE( orthogonality );
    }

    SECTION( "test parity - carthesian coordinates" ) {

        Vector moment1 = testBase.ComputeSphericalBasis( 0, 0 );
        Vector moment2 = testBase.ComputeSphericalBasis( 0, 0 );

        // Parity in carthesian coordinates
        QGaussLegendreTensorized quad( config );

        double x, y, z;

        bool errorWithinBounds = true;
        for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
            x       = quad.GetPoints()[idx_quad][0];
            y       = quad.GetPoints()[idx_quad][1];
            z       = quad.GetPoints()[idx_quad][2];
            moment1 = testBase.ComputeSphericalBasis( x, y, z );
            moment2 = testBase.ComputeSphericalBasis( -x, -y, -z );

            int idx_sys;
            double result = 0.;

            for( int l_idx = 0; l_idx <= int( maxMomentDegree ); l_idx++ ) {
                for( int k_idx = -l_idx; k_idx <= l_idx; k_idx++ ) {
                    idx_sys = testBase.GlobalIdxBasis( l_idx, k_idx );

                    if( l_idx % 2 == 0 )
                        result = moment2[idx_sys] - moment1[idx_sys];
                    else
                        result = moment2[idx_sys] + moment1[idx_sys];

                    if( std::fabs( result ) > 1e2 * std::numeric_limits<double>::epsilon() ) errorWithinBounds = false;
                }
            }
        }
        REQUIRE( errorWithinBounds );
    }

    SECTION( "test parity - polar coordinates" ) {

        Vector moment1 = testBase.ComputeSphericalBasis( 0, 0 );
        Vector moment2 = testBase.ComputeSphericalBasis( 0, 0 );

        int idx_sys;
        double result = 0.;

        // // test in polar coordinates
        bool errorWithinBounds = true;
        for( double my = -1.0; my < 1.0; my += 0.1 ) {
            for( double phi = 0.0; phi < 2 * M_PI; phi += 0.1 ) {
                moment2 = testBase.ComputeSphericalBasis( my, phi );
                moment1 = testBase.ComputeSphericalBasis( -my, M_PI + phi );

                for( int l_idx = 0; l_idx <= int( maxMomentDegree ); l_idx++ ) {
                    for( int k_idx = -l_idx; k_idx <= l_idx; k_idx++ ) {
                        idx_sys = testBase.GlobalIdxBasis( l_idx, k_idx );

                        if( l_idx % 2 == 0 )
                            result = moment2[idx_sys] - moment1[idx_sys];
                        else
                            result = moment2[idx_sys] + moment1[idx_sys];

                        if( std::fabs( result ) > 1e2 * std::numeric_limits<double>::epsilon() ) errorWithinBounds = false;
                    }
                }
            }
        }
        REQUIRE( errorWithinBounds );
    }
}
