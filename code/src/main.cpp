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

int main( int argc, char** argv ) {
    MPI_Init( &argc, &argv );

    // QGaussLegendreTensorized quad( 4 );
    // int maxMomentDegree = 2;
    // SphericalHarmonics testBase( maxMomentDegree );
    //
    // double x, y, z, w;
    // Vector moment = testBase.ComputeSphericalBasis( 0, 1, 0 );
    //// 9 basis moments if degree = 2
    //// unsigned nTotalMoments = moment.size();
    //
    // Vector results( moment.size(), 0.0 );
    //
    // int idx_basis = 0, idx_basisM = 0;
    //
    // for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
    //    x      = quad.GetPoints()[idx_quad][0];
    //    y      = quad.GetPoints()[idx_quad][1];
    //    z      = quad.GetPoints()[idx_quad][2];
    //    w      = quad.GetWeights()[idx_quad];
    //    moment = testBase.ComputeSphericalBasis( x, y, z );
    //
    //    for( int idx_l = 0; idx_l <= 2; idx_l++ ) {
    //        idx_basis = testBase.GlobalIdxBasis( idx_l, 0 );
    //
    //        results[idx_basis] += w * moment[idx_basis] * moment[idx_basis];
    //
    //        for( int idx_k = 1; idx_k <= idx_l; idx_k++ ) {
    //            idx_basis  = testBase.GlobalIdxBasis( idx_l, idx_k );
    //            idx_basisM = testBase.GlobalIdxBasis( idx_l, -idx_k );
    //            // std::cout << "(" << idx_basis << " | " << idx_basisM << ")\n";
    //            results[idx_basis] += w * moment[idx_basis] * moment[idx_basis];
    //            results[idx_basisM] += w * moment[idx_basisM] * moment[idx_basisM];
    //        }
    //    }
    //}
    // std::cout << "Normalization integration:\n " << results << "\n";
    //
    //// Normalized integration
    //// Vector resultsNormal( moment.size(), 0.0 );
    ////
    //// for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
    ////     x      = quad.GetPoints()[idx_quad][0];
    ////     y      = quad.GetPoints()[idx_quad][1];
    ////     z      = quad.GetPoints()[idx_quad][2];
    ////     w      = quad.GetWeights()[idx_quad];
    ////     moment = testBase.ComputeSphericalBasis( x, y, z );
    ////
    ////     for( unsigned idx_sys = 0; idx_sys < 9; idx_sys++ ) {
    ////
    ////         resultsNormal[idx_sys] += w * moment[idx_sys] * moment[idx_sys] / results[idx_sys];
    ////         // std::cout << idx_quad << ": " << results[0] << "\n";
    ////     }
    //// }
    //// std::cout << "moment integration:\n " << resultsNormal << "\n";
    //
    //// Vector results2( moment.size(), 0.0 );
    ////
    //// for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
    ////     x      = quad.GetPoints()[idx_quad][0];
    ////     y      = quad.GetPoints()[idx_quad][1];
    ////     z      = quad.GetPoints()[idx_quad][2];
    ////     w      = quad.GetWeights()[idx_quad];
    ////     moment = testBase.ComputeSphericalBasis( x, y, z );
    ////
    ////     for( unsigned idx_sys = 1; idx_sys < 9; idx_sys++ ) {
    ////         results2[idx_sys] += w * moment[idx_sys - 1] * moment[idx_sys];
    ////         // std::cout << idx_quad << ": " << results[0] << "\n";
    ////     }
    //// }
    //// std::cout << "moment integration:\n " << results2 << "\n";
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

    // Print input file and run info to file
    PrintLogHeader( filename );

    // Build solver
    Solver* solver = Solver::Create( config );

    // Run solver and export
    solver->Solve();
    solver->Save();

    MPI_Finalize();
    return EXIT_SUCCESS;
}
