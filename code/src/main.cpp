#include <Python.h>
#include <mpi.h>

#include "io.h"
#include "solvers/solverbase.h"

#include "settings/config.h"

#include "solvers/sphericalharmonics.h"
#include <fstream>
#include <iostream>
#include <string>

// ----
#include "quadratures/qgausslegendretensorized.h"
#include "quadratures/qmontecarlo.h"
#include "solvers/sphericalharmonics.h"

double testFunc( double my, double phi ) { return my * my + phi; }

double testFunc2( double x, double y, double z ) { return x + y + z; }

int main( int argc, char** argv ) {
    MPI_Init( &argc, &argv );
    wchar_t* program = Py_DecodeLocale( argv[0], NULL );
    Py_SetProgramName( program );

    // QGaussLegendreTensorized quad( 8 );
    // int maxMomentDegree = 8;
    //
    //// test quadrature;
    // std::cout << " Integration test:" << quad.Integrate( testFunc2 ) << "\n";
    //
    // SphericalHarmonics testBase( maxMomentDegree );
    //
    // double my, phi, w, resTest = 0.0;
    // Vector moment = testBase.ComputeSphericalBasis( 0, 1, 0 );
    //// 9 basis moments if degree = 2
    //// unsigned nTotalMoments = moment.size();
    ////
    // Vector results( moment.size(), 0.0 );
    //
    // int idx_basis = 0, idx_basisM = 0;
    //
    // for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
    //    my     = quad.GetPointsSphere()[idx_quad][0];
    //    phi    = quad.GetPointsSphere()[idx_quad][1];
    //    w      = quad.GetWeights()[idx_quad];
    //    moment = testBase.ComputeSphericalBasis( my, phi );
    //
    //    double funcEval = testFunc( my, phi );
    //    std::cout << "f(my= " << my << ", phi= " << phi / M_PI << "*pi) =  " << funcEval << " w " << w << "\n";
    //    resTest += w * funcEval;
    //
    //    // std::cout << "my, phi " << my << " , " << phi << " ::\t mom5 " << moment[5] << " mom7 \t" << moment[7] << "\n";
    //    // results[5] += w * moment[5] * moment[5];
    //    // results[7] += w * moment[7] * moment[7];
    //    //
    //    for( int idx_l = 0; idx_l <= maxMomentDegree; idx_l++ ) {
    //        idx_basis = testBase.GlobalIdxBasis( idx_l, 0 );
    //
    //        results[idx_basis] += w * moment[idx_basis] * moment[idx_basis];
    //
    //        for( int idx_k = 1; idx_k <= idx_l; idx_k++ ) {
    //            idx_basis  = testBase.GlobalIdxBasis( idx_l, idx_k );
    //            idx_basisM = testBase.GlobalIdxBasis( idx_l, -idx_k );
    //
    //            results[idx_basis] += w * moment[idx_basis] * moment[idx_basis];
    //            results[idx_basisM] += w * moment[idx_basisM] * moment[idx_basisM];
    //        }
    //    }
    //}
    // std::cout << " testInt " << resTest << "\n";
    // std::cout << " Squared Integration:\n " << results << "\n";
    //
    // Normalized integration
    // Vector resultsNormal( moment.size(), 0.0 );
    //
    // for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
    //     x      = quad.GetPoints()[idx_quad][0];
    //     y      = quad.GetPoints()[idx_quad][1];
    //     z      = quad.GetPoints()[idx_quad][2];
    //     w      = quad.GetWeights()[idx_quad];
    //     moment = testBase.ComputeSphericalBasis( x, y, z );
    //
    //     for( unsigned idx_sys = 0; idx_sys < 9; idx_sys++ ) {
    //
    //         resultsNormal[idx_sys] += w * moment[idx_sys] * moment[idx_sys] / results[idx_sys];
    //         // std::cout << idx_quad << ": " << results[0] << "\n";
    //     }
    // }
    // std::cout << "moment integration:\n " << resultsNormal << "\n";

    // Vector results2( moment.size(), 0.0 );
    //
    // for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
    //     x      = quad.GetPoints()[idx_quad][0];
    //     y      = quad.GetPoints()[idx_quad][1];
    //     z      = quad.GetPoints()[idx_quad][2];
    //     w      = quad.GetWeights()[idx_quad];
    //     moment = testBase.ComputeSphericalBasis( x, y, z );
    //
    //     for( unsigned idx_sys = 1; idx_sys < 9; idx_sys++ ) {
    //         results2[idx_sys] += w * moment[idx_sys - 1] * moment[idx_sys];
    //         // std::cout << idx_quad << ": " << results[0] << "\n";
    //     }
    // }
    // std::cout << "moment integration:\n " << results2 << "\n";
    //
    // std::ofstream myfile;
    // myfile.open( "assLegendre.csv" );
    //
    // std::vector<double> results3;
    // for( double dMy = -1.0; dMy <= 1.0; dMy += 0.01 ) {
    //    results3 = testBase.GetAssLegendrePoly( dMy );
    //    myfile << dMy;
    //    //     std::cout << dMy;
    //    for( unsigned idx_basis = 0; idx_basis < results3.size(); idx_basis++ ) {
    //        myfile << "," << results3[idx_basis];
    //        // std::cout << " | " << results3[idx_basis];
    //    }
    //    //     std::cout << "\n";
    //    myfile << "\n";
    //}
    // myfile.close();
    //

    std::string filename = ParseArguments( argc, argv );

    // CD  Load Settings from File
    Config* config = new Config( filename );

    // Test the physics reader
    // Physics testPhysic();
    // testPhysic.ReadENDL( "ENDL_H.txt" );

    // Print input file and run info to file
    PrintLogHeader( filename );

    // Build solver
    Solver* solver = Solver::Create( config );

    // Run solver and export
    solver->Solve();
    solver->Save();

    if( Py_IsInitialized() ) Py_Finalize();
    MPI_Finalize();
    return EXIT_SUCCESS;
}
