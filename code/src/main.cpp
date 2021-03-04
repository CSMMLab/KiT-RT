/*! @file: main.cpp
 *  @brief Main method to call the KiT-RT solver suite
 *  @author: J. Kusch, S. Schotthöfer, P. Stammer, J. Wolters, T. Xiao
 *  @version: 0.1
 */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL KITRT_ARRAY_API
#include <mpi.h>
#include <random>
#include <string>

#include "common/config.h"
#include "common/io.h"
#include "solvers/solverbase.h"

#include "toolboxes/datageneratorbase.h"

#ifdef BUILD_GUI
#include <QApplication>

#include "mainwindow.h"
#endif

int main( int argc, char** argv ) {
#ifdef BUILD_GUI
    MPI_Init( &argc, &argv );
    QApplication app( argc, argv );
    MainWindow mw;
    mw.show();
    return app.exec();
#else
    MPI_Init( &argc, &argv );
    wchar_t* program = Py_DecodeLocale( argv[0], NULL );
    Py_SetProgramName( program );

    std::string filename = ParseArguments( argc, argv );

    // CD  Load Settings from File
    Config* config = new Config( filename );

    // Print input file and run info to file
    PrintLogHeader( filename );

    if( config->GetDataGeneratorMode() ) {
        // Build Data generator
        DataGeneratorBase* datagen = DataGeneratorBase::Create( config );
        // Generate Data and export
        datagen->ComputeTrainingData();
    }
    else {
        // Build solver
        unsigned nMCSamples = 10;
        std::vector<double> expDosis;
        std::random_device rd;       // Will be used to obtain a seed for the random number engine
        std::mt19937 gen( rd() );    // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis( 1.0, 2.0 );
        for( unsigned k = 0; k < nMCSamples; ++k ) {
            SolverBase* solver = SolverBase::Create( config );
            double density     = dis( gen );
            solver->SetDensity( density );    // if density < 1, then you have to adapt the CFL condition, use density \in [1,2]

            // Run solver and export
            solver->Solve();
            std::vector<double> dosis = solver->GetDosis();

            // update expected dosis
            if( k == 0 ) expDosis = std::vector<double>( solver->GetNCells(), 0.0 );
            for( unsigned j = 0; j < config->GetNCells(); ++j ) {
                expDosis[j] = expDosis[j] + dosis[j] / nMCSamples;
            }
            delete solver;
        }
        // write expected dosis as vtk file
        std::vector<std::string> fieldNames{ "dose", "normalized dose" };
        std::vector<std::vector<std::string>> fieldNamesWrapper{ fieldNames };
        std::vector<std::vector<double>> dose( 1, expDosis );
        std::vector<std::vector<double>> normalizedDose( 1, expDosis );
        double maxDose = *std::max_element( expDosis.begin(), expDosis.end() );
        for( unsigned i = 0; i < expDosis.size(); ++i ) normalizedDose[0][i] /= maxDose;
        std::vector<std::vector<std::vector<double>>> results{ dose, normalizedDose };
        ExportVTK( config->GetOutputFile(), results, fieldNamesWrapper, LoadSU2MeshFromFile( config ) );
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
#endif
}
