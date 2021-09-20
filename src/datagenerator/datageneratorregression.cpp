#include "datagenerator/datageneratorregression.h"
#include "common/config.h"
#include "entropies/entropybase.h"
#include "optimizers/newtonoptimizer.h"
#include "toolboxes/errormessages.h"

#include <math.h>

DataGeneratorRegression::DataGeneratorRegression( Config* settings ) : DataGeneratorBase( settings ) {
    _gridSize = _setSize;    // Number of Datapoints in first dimension
}

DataGeneratorRegression::~DataGeneratorRegression() {}

void DataGeneratorRegression::ComputeTrainingData() {
    PrintLoadScreen();

    auto log = spdlog::get( "event" );
    if( _settings->GetAlphaSampling() ) {
        // --- sample alpha ---
        log->info( "| Sample Lagrange multipliers." );
        SampleMultiplierAlpha();
        log->info( "| Multipliers sampled." );

        log->info( "| Compute realizable problems." );

        // --- Postprocessing
        ComputeRealizableSolution();
    }
    else {
        // --- sample u ---
        SampleSolutionU();
        log->info( "| Moments sampled." );
        log->info( "| Start solving the optimization problems. This may take some minutes." );

        // ---- Check realizability ---
        // CheckRealizability();
        /*
        std::vector<double> timings( 100, 0.0 );

        for( unsigned j = 0; j < _setSize; j++ ) {
            _uSol[j][0] = 1.0;
            _uSol[j][1] = 0.99;
            _uSol[j][2] = 0.99 * 0.99;
        }

         for( unsigned i = 0; i < 100; i++ ) {
        // Record start time
        auto start = std::chrono::high_resolution_clock::now();
        // --- compute alphas ---
        _optimizer->SolveMultiCell( _alpha, _uSol, _momentBasis );

        // Record end time
        auto finish                           = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = finish - start;
        timings[i]                            = elapsed.count();
        log->info( "| Elapsed time for solving the entropy minimization problem: " + std::to_string( elapsed.count() ) + " seconds. Iteration" +
                   std::to_string( i ) + "." );
        // reset training solution
        for( unsigned j = 0; j < _setSize; j++ ) {
            for( unsigned k = 0; k < _nTotalEntries; k++ ) {
                _alpha[j][k] = 0.0;
            }
        }
        }

        // Compute statistics
        double mean = 0;
        for( unsigned i = 0; i < 100; i++ ) {
            mean += timings[i];
        }
        mean /= 100;
        double stdDev = 0;
        for( unsigned i = 0; i < 100; i++ ) {
            stdDev += ( timings[i] - mean ) * ( timings[i] - mean );
        }
        stdDev /= 100;
        stdDev = sqrt( stdDev );
        log->info( "| Mean timing: " + std::to_string( mean ) + " seconds. Standard deviation: " + std::to_string( stdDev ) + " seconds." );
        */

        // some checks
        // Vector u     = { 1.0, 0.0, 0.5, 0.0 };
        // Vector alpha = { 1.0, 0.0, 0.5, 0.01 };
        //_optimizer->Solve( alpha, u, _momentBasis, 0 );
        // std::cout << "u=" << u << "\nalpha=" << alpha << "\n";
        // u     = { 1.0, 0.0, 0.5, 0.01 };
        // alpha = { 1.0, 0.0, 0.5, 0.01 };
        //_optimizer->Solve( alpha, u, _momentBasis, 0 );
        // std::cout << "u=" << u << "\nalpha=" << alpha << "\n";
        // u     = { 1.0, 0.0, 0.5, 0.05 };
        // alpha = { 1.0, 0.0, 0.5, 0.01 };
        //_optimizer->Solve( alpha, u, _momentBasis, 0 );
        // std::cout << "u=" << u << "\nalpha=" << alpha << "\n";
        // u     = { 1.0, 0.0, 0.5, 0.1 };
        // alpha = { 1.0, 0.0, 0.5, 0.01 };
        //_optimizer->Solve( alpha, u, _momentBasis, 0 );
        // std::cout << "u=" << u << "\nalpha=" << alpha << "\n";
        // u     = { 1.0, 0.0, 0.5, 0.2 };
        // alpha = { 1.0, 0.0, 0.5, 0.01 };
        //_optimizer->Solve( alpha, u, _momentBasis, 0 );
        // std::cout << "u=" << u << "\nalpha=" << alpha << "\n";

        // --- compute alphas ---
        _optimizer->SolveMultiCell( _alpha, _uSol, _momentBasis );
        // --- Postprocessing
        if( _settings->GetRelizabilityReconsU() ) {
            log->info( "| Compute realizable problems." );
            ComputeRealizableSolution();
        }
    }

    log->info( "| Compute entropies." );

    // --- compute entropy functional ---
    ComputeEntropyH_primal();

    log->info( "| Print Solution." );

    // --- Print everything ----
    PrintTrainingData();
}

void DataGeneratorRegression::ComputeEntropyH_dual() {
#pragma omp parallel for schedule( guided )
    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        _hEntropy[idx_set] = _optimizer->ComputeObjFunc( _alpha[idx_set], _uSol[idx_set], _momentBasis );
    }
}

void DataGeneratorRegression::ComputeEntropyH_primal() {
#pragma omp parallel for schedule( guided )
    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        double result = 0.0;
        // Integrate (eta(eta'_*(alpha*m))
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            result += _entropy->Entropy( _entropy->EntropyPrimeDual( dot( _alpha[idx_set], _momentBasis[idx_quad] ) ) ) * _weights[idx_quad];
        }
        _hEntropy[idx_set] = result;
    }

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        if( isnan( _hEntropy[idx_set] ) ) {
            std::string msgU     = "(";
            std::string msgAlpha = "(";
            for( unsigned idx_basis = 0; idx_basis < _nTotalEntries - 1; idx_basis++ ) {
                msgU += std::to_string( _uSol[idx_set][idx_basis] ) + "|";
                msgAlpha += std::to_string( _alpha[idx_set][idx_basis] ) + "|";
            }
            msgU += std::to_string( _uSol[idx_set][_nTotalEntries - 1] ) + ")";
            msgAlpha += std::to_string( _uSol[idx_set][_nTotalEntries - 1] ) + ")";

            ErrorMessages::Error( "Value for h is NaN. This can happen near the boundary of the realizable set.\nu= " + msgU +
                                      "\nalpha= " + msgAlpha +
                                      "\nPlease adjust the Options "
                                      "REALIZABLE_SET_EPSILON_U0 and REALIZABLE_SET_EPSILON_U1.",
                                  CURRENT_FUNCTION );
        }
    }
}

void DataGeneratorRegression::PrintTrainingData() {
    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );
    log->info( "---------------------- Data Generation Successful ------------------------" );

    std::string uSolString  = "";
    std::string alphaString = "";
    for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
        uSolString += "u_" + std::to_string( idx_sys ) + ",";
        alphaString += "alpha_" + std::to_string( idx_sys ) + ",";
    }
    // log->info( uSolString + alphaString + "h" );
    logCSV->info( uSolString + alphaString + "h" );

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {

        std::stringstream streamU, streamAlpha, streamH;

        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            streamU << std::fixed << std::setprecision( 12 ) << _uSol[idx_set][idx_sys] << ",";
            streamAlpha << std::fixed << std::setprecision( 12 ) << _alpha[idx_set][idx_sys] << ",";
        }
        streamH << std::fixed << std::setprecision( 12 ) << _hEntropy[idx_set];

        std::string uSolString  = streamU.str();
        std::string alphaString = streamAlpha.str();
        std::string hString     = streamH.str();

        // log->info(  uSolString + alphaString + hString  );
        logCSV->info( uSolString + alphaString + hString );
    }
}

void DataGeneratorRegression::ComputeSetSizeAlpha() {
    if( _settings->GetSizeByDimension() ) {
        _setSize = _gridSize;

        for( unsigned i = 0; i < _nTotalEntries - 2; i++ ) {
            _setSize *= _gridSize;
        }
    }
}
