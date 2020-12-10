/*!
 * \file datagenerator.cpp
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#include "toolboxes/datagenerator.h"
#include "common/config.h"
#include "entropies/entropybase.h"
#include "optimizers/newtonoptimizer.h"
#include "quadratures/quadraturebase.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/sphericalbase.h"

#include "spdlog/spdlog.h"

#include <iostream>

nnDataGenerator::nnDataGenerator( Config* settings ) {
    _settings = settings;
    _setSize  = settings->GetTrainingDataSetSize();

    _LMaxDegree = settings->GetMaxMomentDegree();

    // Quadrature
    _quadrature       = QuadratureBase::Create( settings );
    _nq               = _quadrature->GetNq();
    _quadPoints       = _quadrature->GetPoints();
    _weights          = _quadrature->GetWeights();
    _quadPointsSphere = _quadrature->GetPointsSphere();
    _settings->SetNQuadPoints( _nq );

    // Spherical Harmonics
    if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS && _LMaxDegree > 0 ) {
        ErrorMessages::Error( "No sampling algorithm for spherical harmonics basis with degree higher than 0 implemented", CURRENT_FUNCTION );
    }
    _basis = SphericalBase::Create( _settings );

    _nTotalEntries = _basis->GetBasisSize();

    _moments = VectorVector( _nq, Vector( _nTotalEntries, 0.0 ) );

    ComputeMoments();

    // Optimizer
    _optimizer = new NewtonOptimizer( _settings );

    // Entropy
    _entropy = EntropyBase::Create( _settings );

    // Initialize Training Data
    _uSol     = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _alpha    = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _hEntropy = std::vector<double>( _setSize, 0.0 );
}

void nnDataGenerator::computeTrainingData() {
    // Prototype: Only for _LMaxDegree == 1
    // Prototype: u is sampled from [0,100]

    // --- sample u ---
    sampleSolutionU();

    // --- compute alphas ---
    _optimizer->SolveMultiCell( _alpha, _uSol, _moments );

    // --- compute entropy functional ---
    computeEntropyH_primal();

    // --- Print everything ----
    printTrainingData();
}

int nnDataGenerator::GlobalIndex( int l, int k ) const {
    int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
    int prevIndicesThisLevel = k + l;    // number of previous indices in current level
    return numIndicesPrevLevel + prevIndicesThisLevel;
}

void nnDataGenerator::ComputeMoments() {
    double my, phi;

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my  = _quadPointsSphere[idx_quad][0];
        phi = _quadPointsSphere[idx_quad][1];

        _moments[idx_quad] = _basis->ComputeSphericalBasis( my, phi );
    }
}

void nnDataGenerator::sampleSolutionU() {
    double du = 100.0 / (double)_setSize;    // Prototype: u is sampled from [0,100]

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        _uSol[idx_set][0] = du * idx_set;
    }
}

void nnDataGenerator::computeEntropyH_dual() {
    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        _hEntropy[idx_set] = _optimizer->ComputeObjFunc( _alpha[idx_set], _uSol[idx_set], _moments );
    }
}

void nnDataGenerator::computeEntropyH_primal() {
    double result = 0.0;

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        result = 0.0;
        // Integrate (eta(eta'_*(alpha*m))
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            result += _entropy->Entropy( _entropy->EntropyPrimeDual( dot( _alpha[idx_set], _moments[idx_quad] ) ) ) * _weights[idx_quad];
        }
        _hEntropy[idx_set] = result;
    }
}

void nnDataGenerator::printTrainingData() {
    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        log->info( "{}, {}, {}", _uSol[idx_set][0], _alpha[idx_set][0], _hEntropy[idx_set] );
        logCSV->info( "{}, {}, {}", _uSol[idx_set][0], _alpha[idx_set][0], _hEntropy[idx_set] );
    }
}
