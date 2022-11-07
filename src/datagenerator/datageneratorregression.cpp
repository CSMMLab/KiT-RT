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
    }
    else {
        // --- sample u ---
        SampleSolutionU();
        log->info( "| Moments sampled." );
        log->info( "| Start solving the optimization problems. This may take some minutes." );

        // ---- Check realizability ---
        // CheckRealizability();

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
