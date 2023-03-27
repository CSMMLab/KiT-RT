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
        Vector alpha_norm_dummy( _setSize, 0 );

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
        _optimizer->SolveMultiCell( rot_alpha_comp, rot_uSol, _momentBasis, alpha_norm_dummy );

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

        // --- compute alphas ---
        Vector alpha_norm_dummy( _setSize, 0 );

        _optimizer->SolveMultiCell( _alpha, _uSol, _momentBasis, alpha_norm_dummy );
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
    if( _reducedSampling ) {
        VectorVector momentsRed = VectorVector( _nq, Vector( _nTotalEntries - 1, 0.0 ) );

        for( unsigned idx_nq = 0; idx_nq < _nq; idx_nq++ ) {    // copy (reduced) moments
            for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                momentsRed[idx_nq][idx_sys - 1] = _momentBasis[idx_nq][idx_sys];
            }
        }
#pragma omp parallel for schedule( guided )
        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            Vector uSolRed( _nTotalEntries - 1, 0.0 );
            Vector alphaRed( _nTotalEntries - 1, 0.0 );
            for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                alphaRed[idx_sys - 1] = _alpha[idx_set][idx_sys];
                uSolRed[idx_sys - 1]  = _uSol[idx_set][idx_sys];
            }
            _hEntropy[idx_set] = -1 * _optimizer->ComputeObjFunc( alphaRed, uSolRed, momentsRed );
        }
    }
    else {
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
