/*!
 * \file datageneratorbase.cpp
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#include "toolboxes/datageneratorbase.h"
#include "common/config.h"
#include "entropies/entropybase.h"
#include "optimizers/newtonoptimizer.h"
#include "quadratures/quadraturebase.h"
#include "spdlog/spdlog.h"
#include "toolboxes/datagenerator1D.h"
#include "toolboxes/datagenerator2D.h"
#include "toolboxes/datagenerator3D.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/sphericalbase.h"
#include "toolboxes/textprocessingtoolbox.h"

#include <iomanip>
#include <math.h>
#include <omp.h>
#include <sstream>

DataGeneratorBase::DataGeneratorBase( Config* settings ) {
    _settings = settings;
    _setSize  = settings->GetTrainingDataSetSize();
    _gridSize = _setSize;

    _maxPolyDegree = settings->GetMaxMomentDegree();

    // Check consistency between dimension of quadrature and sample basis
    if( _settings->GetDim() == 1 ) {
        if( _settings->GetQuadName() != QUAD_GaussLegendre1D && _settings->GetQuadName() != QUAD_GaussChebyshev1D ) {
            ErrorMessages::Error( "For 1D Sampling, please choose a 1D quadrature rule.", CURRENT_FUNCTION );
        }
    }
    else {
        if( _settings->GetQuadName() == QUAD_GaussLegendre1D || _settings->GetQuadName() == QUAD_GaussChebyshev1D ) {
            ErrorMessages::Error( "For 3D Sampling, please choose a 3D quadrature rule.", CURRENT_FUNCTION );
        }
    }
    // Quadrature
    _quadrature       = QuadratureBase::Create( settings );
    _nq               = _quadrature->GetNq();
    _quadPoints       = _quadrature->GetPoints();
    _weights          = _quadrature->GetWeights();
    _quadPointsSphere = _quadrature->GetPointsSphere();

    // Spherical Harmonics
    if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS && _maxPolyDegree > 0 ) {
        ErrorMessages::Error( "No sampling algorithm for spherical harmonics basis with degree higher than 0 implemented", CURRENT_FUNCTION );
    }
    _basis = SphericalBase::Create( _settings );

    _nTotalEntries = _basis->GetBasisSize();

    _moments = VectorVector( _nq, Vector( _nTotalEntries, 0.0 ) );

    // Optimizer
    _optimizer = new NewtonOptimizer( _settings );

    // Entropy
    _entropy = EntropyBase::Create( _settings );
}

DataGeneratorBase::~DataGeneratorBase() {
    delete _quadrature;
    delete _entropy;
}

DataGeneratorBase* DataGeneratorBase::Create( Config* settings ) {
    switch( settings->GetDim() ) {
        case 1: return new DataGenerator1D( settings );
        case 2: return new DataGenerator2D( settings );
        case 3: return new DataGenerator3D( settings );
        default: ErrorMessages::Error( "Sampling for more than 3 dimensions is not yet supported.", CURRENT_FUNCTION );
    }
    return nullptr;
}

void DataGeneratorBase::ComputeTrainingData() {
    PrintLoadScreen();
    auto log = spdlog::get( "event" );

    if( _settings->GetAlphaSampling() ) {
        // --- sample alpha ---
        SampleMultiplierAlpha();
        log->info( "| Multipliers sampled." );

        log->info( "| Making moments realizable problems." );

        // --- Postprocessing
        ComputeRealizableSolution();
    }
    else {
        // --- sample u ---
        SampleSolutionU();
        log->info( "| Moments sampled." );
        log->info( "| Start solving the optimization problems. This may take some minutes." );

        // ---- Check realizability ---
        CheckRealizability();

        // --- compute alphas ---

        _optimizer->SolveMultiCell( _alpha, _uSol, _moments );

        log->info( "| Making moments realizable problems." );

        // --- Postprocessing
        if( _settings->GetRelizabilityReconsU() ) {
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

void DataGeneratorBase::SampleMultiplierAlpha() {
    double minAlphaValue = -50;
    double maxAlphaValue = 50;

    if( _settings->GetNormalizedSampling() ) {
        // compute reduced version of alpha and m
        if( _maxPolyDegree == 0 ) {
            ErrorMessages::Error( "Normalized sampling not meaningful for M0 closure", CURRENT_FUNCTION );
        }
        VectorVector alphaRed   = VectorVector( _setSize, Vector( _nTotalEntries - 1, 0.0 ) );
        VectorVector momentsRed = VectorVector( _nq, Vector( _nTotalEntries - 1, 0.0 ) );

        for( unsigned idx_nq = 0; idx_nq < _nq; idx_nq++ ) {    // copy (reduced) moments
            for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                momentsRed[idx_nq][idx_sys - 1] = _moments[idx_nq][idx_sys];
            }
        }
        double dalpha = 0;
        switch( _maxPolyDegree ) {
            case 1:
                // Sample alpha1 from [minAlphaValue, maxAlphaValue], then compute alpha_0 s.t. u_0 = 1
                dalpha = ( maxAlphaValue - minAlphaValue ) / (double)_setSize;

                for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
                    alphaRed[idx_set][0] = minAlphaValue + idx_set * dalpha;
                }
                break;
            default: ErrorMessages::Error( "Not yet implemented!", CURRENT_FUNCTION );
        }

        // Compute alpha_0 = log(<exp(alpha m )>) // for maxwell boltzmann! only
        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            double integral = 0.0;
            // Integrate (eta(eta'_*(alpha*m))
            for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                integral += _entropy->EntropyPrimeDual( dot( alphaRed[idx_set], momentsRed[idx_quad] ) ) * _weights[idx_quad];
            }
            _alpha[idx_set][0] = -log( integral );    // log trafo

            // copy all other alphas to the member
            for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                _alpha[idx_set][idx_sys] = alphaRed[idx_set][idx_sys - 1];
            }
        }
    }
    else {
        // non normalized sampling
        // TODO
        ErrorMessages::Error( "Not yet implemented!", CURRENT_FUNCTION );
    }
}

void DataGeneratorBase::ComputeEntropyH_dual() {
#pragma omp parallel for schedule( guided )
    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        _hEntropy[idx_set] = _optimizer->ComputeObjFunc( _alpha[idx_set], _uSol[idx_set], _moments );
    }
}

void DataGeneratorBase::ComputeEntropyH_primal() {
#pragma omp parallel for schedule( guided )
    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        double result = 0.0;
        // Integrate (eta(eta'_*(alpha*m))
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            result += _entropy->Entropy( _entropy->EntropyPrimeDual( dot( _alpha[idx_set], _moments[idx_quad] ) ) ) * _weights[idx_quad];
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

void DataGeneratorBase::PrintTrainingData() {
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

void DataGeneratorBase::ComputeRealizableSolution() {
#pragma omp parallel for schedule( guided )
    for( unsigned idx_sol = 0; idx_sol < _setSize; idx_sol++ ) {
        double entropyReconstruction = 0.0;
        _uSol[idx_sol]               = 0;
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            // Make entropyReconstruction a member vector, s.t. it does not have to be re-evaluated in ConstructFlux
            entropyReconstruction = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_sol], _moments[idx_quad] ) );
            _uSol[idx_sol] += _moments[idx_quad] * ( _weights[idx_quad] * entropyReconstruction );
        }
    }
}

void DataGeneratorBase::PrintLoadScreen() {
    auto log = spdlog::get( "event" );
    log->info( "------------------------ Data Generation Starts --------------------------" );
    log->info( "| Generating {} datapoints.", _setSize );
}

void DataGeneratorBase::AdaptBasisSize() {
    // Remove zero order Moment for dimension reduction

    if( _settings->GetNormalizedSampling() ) {
        VectorVector momentTemp = _moments;
        _nTotalEntries          = _nTotalEntries - 1;
        _moments                = VectorVector( _nq, Vector( _nTotalEntries, 0.0 ) );

        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
                _moments[idx_quad][idx_sys] = momentTemp[idx_quad][idx_sys + 1];
            }
        }
    }
}
