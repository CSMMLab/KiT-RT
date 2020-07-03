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
    //
    // QGaussLegendreTensorized quad( 4 );
    // SphericalHarmonics testBase( 2 );
    //
    // double x, y, z, w;
    // Vector moment = testBase.ComputeSphericalBasis( 0, 1, 0 );
    // // 9 basis moments if degree = 2
    //
    // Vector results( moment.size(), 0.0 );
    //
    // for( unsigned idx_quad = 0; idx_quad < quad.GetNq(); idx_quad++ ) {
    //     x      = quad.GetPoints()[idx_quad][0];
    //     y      = quad.GetPoints()[idx_quad][1];
    //     z      = quad.GetPoints()[idx_quad][2];
    //     w      = quad.GetWeights()[idx_quad];
    //     moment = testBase.ComputeSphericalBasis( x, y, z );
    //
    //     for( unsigned idx_sys = 0; idx_sys < 9; idx_sys++ ) {
    //         results[idx_sys] += w * moment[idx_sys] * moment[idx_sys];
    //         // std::cout << idx_quad << ": " << results[0] << "\n";
    //     }
    // }
    // std::cout << "moment integration:\n " << results << "\n";
    //
    // // Normalized integration
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
    //
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

    std::string filename = ParseArguments( argc, argv );

    // Load Settings from File
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
