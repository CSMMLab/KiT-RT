#include <mpi.h>

#include "io.h"
#include "solvers/solverbase.h"

#include "settings/config.h"

#include "solvers/sphericalharmonics.h"
#include <fstream>
#include <iostream>
#include <string>

int main( int argc, char** argv ) {
    MPI_Init( &argc, &argv );

    //    std::string filename = ParseArguments( argc, argv );
    //
    //    // Load Settings from File
    //    Config* config = new Config( filename );
    //
    //    // Print input file and run info to file
    //    PrintLogHeader( filename );
    //
    //    // Build solver
    //    Solver* solver = Solver::Create( config );
    //
    //    // Run solver and export
    //    solver->Solve();
    //    solver->Save();
    //

    SphericalHarmonics testBase( 2 );
    double d_my  = 2. / 50.;
    double d_phi = 2. * M_PI / 50.;
    std::vector<double> erg( 4, 0.0 );

    std::ofstream outputFile;
    std::string filename = "harmonicTest.csv";
    outputFile.open( filename );

    outputFile << "my"
               << ","
               << "phi"
               << ","
               << "v0"
               << ","
               << "v1"
               << ","
               << "v2"
               << ","
               << "v3"
               << ","
               << "v4" << std::endl;

    for( int my_idx = 0; my_idx < 50; my_idx++ ) {
        for( int phi_idx = 0; phi_idx < 50; phi_idx++ ) {
            erg = testBase.ComputeSphericalBasis( my_idx * d_my, phi_idx * d_phi );
            outputFile << my_idx * d_my << "," << phi_idx * d_phi << "," << erg[0] << "," << erg[1] << "," << erg[2] << "," << erg[3] << "," << erg[4]
                       << std::endl;
        }
    }

    outputFile.close();

    MPI_Finalize();
    return EXIT_SUCCESS;
}
