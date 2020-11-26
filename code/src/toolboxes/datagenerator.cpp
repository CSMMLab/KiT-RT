/*!
 * \file datagenerator.cpp
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#include "toolboxes/datagenerator.h"
#include "common/config.h"
#include "quadratures/quadraturebase.h"
#include "solvers/sphericalharmonics.h"

#include "spdlog/spdlog.h"

#include <iostream>

nnDataGenerator::nnDataGenerator( Config* settings ) {
    _settings = settings;
    _setSize  = settings->GetTrainingDataSetSize();

    _LMaxDegree    = settings->GetMaxMomentDegree();
    _nTotalEntries = (unsigned)GlobalIndex( _LMaxDegree, _LMaxDegree ) + 1;

    // build quadrature object and store frequently used params
    _quadrature       = QuadratureBase::CreateQuadrature( settings );
    _nq               = _quadrature->GetNq();
    _quadPoints       = _quadrature->GetPoints();
    _weights          = _quadrature->GetWeights();
    _quadPointsSphere = _quadrature->GetPointsSphere();
    _settings->SetNQuadPoints( _nq );

    _basis   = new SphericalHarmonics( _LMaxDegree );
    _moments = VectorVector( _nq, Vector( _nTotalEntries, 0.0 ) );
    ComputeMoments();

    // Initialize Training Data
    _uSol     = std::vector( _setSize, std::vector( _nTotalEntries, 0.0 ) );
    _alpha    = std::vector( _setSize, std::vector( _nTotalEntries, 0.0 ) );
    _hEntropy = std::vector( _setSize, 0.0 );
}

void nnDataGenerator::computeTrainingData() {
    // Prototype: Only for _LMaxDegree == 1
    // Prototype: u is sampled from [0,100]

    // --- sample u ---
    sampleSolutionU();
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
    auto log  = spdlog::get( "event" );

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        _uSol[idx_set][0] = du * idx_set;
        log->info( "{}", _uSol[idx_set][0] );
    }
    log->flush();
}
