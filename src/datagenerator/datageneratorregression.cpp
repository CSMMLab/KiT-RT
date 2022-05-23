/*!
 * \file datageneratorregression.h
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#include "datagenerator/datageneratorregression.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "optimizers/newtonoptimizer.hpp"
#include "toolboxes/errormessages.hpp"

#include "toolboxes/textprocessingtoolbox.hpp"

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

        // debugging purposes for rotation
        TextProcessingToolbox::PrintVectorVector( _uSol );
        TextProcessingToolbox::PrintVectorVector( _alpha );

        std::cout << "here\n";

        VectorVector rot_uSol       = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
        VectorVector rot_alpha_comp = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
        VectorVector rot_alpha      = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );

        for( unsigned i = 0; i < _setSize; i++ ) {
            Vector u1{ _uSol[i][1], _uSol[i][2] };
            Matrix u2{ { _uSol[i][3], _uSol[i][4] }, { _uSol[i][4], _uSol[i][5] } };

            Vector alpha1{ _alpha[i][1], _alpha[i][2] };
            Matrix alpha2{ { _alpha[i][3], _alpha[i][4] }, { _alpha[i][4], _alpha[i][5] } };

            Matrix rotationMat  = CreateRotator( u1 );
            Matrix rotationMatT = blaze::trans( rotationMat );

            u1 = RotateM1( u1, rotationMat );
            u2 = RotateM2( u2, rotationMat, rotationMatT );

            rot_uSol[i][0] = (float)( u1[0] );
            rot_uSol[i][1] = (float)( u1[1] );    // should be zero
            rot_uSol[i][2] = (float)( u2( 0, 0 ) );
            rot_uSol[i][3] = (float)( u2( 0, 1 ) );
            rot_uSol[i][4] = (float)( u2( 1, 1 ) );

            rot_alpha[i][0] = (float)( alpha1[0] );
            rot_alpha[i][1] = (float)( alpha1[1] );    // should be zero
            rot_alpha[i][2] = (float)( alpha2( 0, 0 ) );
            rot_alpha[i][3] = (float)( alpha2( 0, 1 ) );
            rot_alpha[i][4] = (float)( alpha2( 1, 1 ) );
        }
        _optimizer->SolveMultiCell( rot_alpha_comp, rot_uSol, _momentBasis );

        TextProcessingToolbox::PrintVectorVector( rot_alpha );
        TextProcessingToolbox::PrintVectorVector( rot_alpha_comp );
        std::cout << "here\n";
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
        if( _settings->GetRealizabilityReconstruction() ) {
            log->info( "| Compute realizable problems." );
            ComputeRealizableSolution();
        }
    }

    log->info( "| Compute entropies." );

    // --- compute entropy functional ---
    // ComputeEntropyH_dual();
    ComputeEntropyH_primal();

    log->info( "| Print Solution." );

    // --- Print everything ----
    PrintTrainingData();
}

Matrix DataGeneratorRegression::CreateRotator( const Vector& uFirstMoment ) {
    double a = uFirstMoment[0];
    double b = uFirstMoment[1];
    double c, s, r;

    r = norm( uFirstMoment );    // sqrt( a * a + b * b );
    c = a / r;
    s = -b / r;

    return Matrix{ { c, -s }, { s, c } };    // Rotation Matrix
}

Vector DataGeneratorRegression::RotateM1( Vector& vec, Matrix& R ) { return R * vec; }

Matrix DataGeneratorRegression::RotateM2( Matrix& vec, Matrix& R, Matrix& Rt ) { return R * vec * Rt; }

void DataGeneratorRegression::ComputeEntropyH_dual() {
#pragma omp parallel for schedule( guided )
    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        _hEntropy[idx_set] = _optimizer->ComputeObjFunc( _alpha[idx_set], _uSol[idx_set], _momentBasis );
    }
}

void DataGeneratorRegression::ComputeEntropyH_primal() {
#pragma omp parallel for schedule( guided )
    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        _hEntropy[idx_set] = -1 * _optimizer->ComputeObjFunc( _alpha[idx_set], _uSol[idx_set], _momentBasis );
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
