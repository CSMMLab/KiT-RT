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

#include <math.h>
#include <omp.h>

DataGeneratorBase::DataGeneratorBase( Config* settings ) {
    _settings = settings;
    _setSize  = settings->GetTrainingDataSetSize();
    _gridSize = _setSize;

    _LMaxDegree = settings->GetMaxMomentDegree();

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
    if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS && _LMaxDegree > 0 ) {
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

    // --- sample u ---
    SampleSolutionU();
    auto log = spdlog::get( "event" );
    log->info( "| Moments sampled." );
    log->info( "| Start solving the optimization problems. This may take some minutes." );

    // ---- Check realizability ---
    CheckRealizability();

    // --- compute alphas ---
    _optimizer->SolveMultiCell( _alpha, _uSol, _moments );

    log->info( "| Making moments realizable problems." );

    // --- Postprocessing
    ComputeRealizableSolution();

    log->info( "| Compute entropies." );

    // --- compute entropy functional ---
    ComputeEntropyH_primal();

    log->info( "| Print Solution." );

    // --- Print everything ----
    PrintTrainingData();
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
        std::string uSolString  = "";
        std::string alphaString = "";
        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            uSolString += std::to_string( _uSol[idx_set][idx_sys] ) + ",";
            alphaString += std::to_string( _alpha[idx_set][idx_sys] ) + ",";
        }
        // log->info( uSolString + alphaString + "{}", _hEntropy[idx_set] );
        logCSV->info( uSolString + alphaString + "{}", _hEntropy[idx_set] );
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
