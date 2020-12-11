/*!
 * \file datagenerator.cpp
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#include "toolboxes/datagenerator.h"
#include "common/config.h"
#include "entropies/entropybase.h"
#include "optimizers/newtonoptimizer.h"
#include "quadratures/qlebedev.h"
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
    SampleSolutionU();

    // --- compute alphas ---
    _optimizer->SolveMultiCell( _alpha, _uSol, _moments );

    // --- compute entropy functional ---
    ComputeEntropyH_primal();

    // --- Print everything ----
    PrintTrainingData();
}

void nnDataGenerator::ComputeMoments() {
    double my, phi;

    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my  = _quadPointsSphere[idx_quad][0];
        phi = _quadPointsSphere[idx_quad][1];

        _moments[idx_quad] = _basis->ComputeSphericalBasis( my, phi );
    }
}

void nnDataGenerator::SampleSolutionU() {
    // Use necessary conditions from Monreal, Dissertation, Chapter 3.2.1, Page 26

    // --- Determine stepsizes etc ---
    double du0 = 10.0 / (double)_setSize;    // Prototype: u is sampled from [0,10]

    // different processes for different
    if( _LMaxDegree == 0 ) {
        // --- sample u in order 0 ---
        // u_0 = <1*psi>

        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            _uSol[idx_set][0] = du0 * idx_set;
        }
    }
    if( _LMaxDegree == 1 ) {
        // Sample points on unit sphere. (Use Lebedev quadrature)
        unsigned order       = 5;    // is avail.
        QLebedev* quad       = new QLebedev( order );
        VectorVector qpoints = quad->GetPoints();    // carthesian coordinates.
        unsigned long nq     = (unsigned long)quad->GetNq();

        // Allocate memory.
        unsigned long setSizeU0 = _setSize;
        unsigned long setSizeU1 = nq * setSizeU0 * ( setSizeU0 + 1 ) / 2;
        _setSize                = setSizeU1;

        // REFACTOR THIS
        _uSol     = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
        _alpha    = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
        _hEntropy = std::vector<double>( _setSize, 0.0 );

        // --- sample u in order 0 ---
        // u_0 = <1*psi>

        for( unsigned long idx_set = 0; idx_set < setSizeU0; idx_set++ ) {
            unsigned long outerIdx = ( idx_set - 1 ) * ( idx_set ) / 2;    // sum over all radii up to current
            outerIdx *= nq;                                                // in each radius step, use all quad points

            //  --- sample u in order 1 ---
            /* order 1 has 3 elements. (omega_x, omega_y,  omega_z) = omega, let u_1 = (u_x, u_y, u_z) = <omega*psi>
             * Condition u_0 >= norm(u_1)   */

            // loop over all radii
            for( unsigned long idx_subset = 0; idx_subset < idx_set; idx_subset++ ) {
                double radius          = du0 * ( idx_subset );    // dont use radius 0
                unsigned long localIdx = outerIdx + idx_subset * nq;

                for( unsigned quad_idx = 0; quad_idx < nq; quad_idx++ ) {
                    unsigned long radiusIdx = localIdx + quad_idx;    // gives the global index

                    _uSol[radiusIdx][0] = radius;
                    // scale quadpoints with radius
                    _uSol[radiusIdx][1] = radius * qpoints[quad_idx][0];
                    _uSol[radiusIdx][2] = radius * qpoints[quad_idx][1];
                    _uSol[radiusIdx][3] = radius * qpoints[quad_idx][2];
                    std::cout << " radiusIdx: " << radiusIdx << "\n";
                    // std::cout << " _uSol[radiusIdx][0]: " << _uSol[radiusIdx][0] << "\n";
                    // std::cout << " _uSol[radiusIdx][1]: " << _uSol[radiusIdx][1] << "\n";
                    // std::cout << " _uSol[radiusIdx][2]: " << _uSol[radiusIdx][2] << "\n";
                    // std::cout << " _uSol[radiusIdx][3]: " << _uSol[radiusIdx][3] << "\n";
                }
            }
        }
    }
    else {
        ErrorMessages::Error( "Sampling for order higher than 1 is not yet supported", CURRENT_FUNCTION );
    }
}

void nnDataGenerator::ComputeEntropyH_dual() {
    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        _hEntropy[idx_set] = _optimizer->ComputeObjFunc( _alpha[idx_set], _uSol[idx_set], _moments );
    }
}

void nnDataGenerator::ComputeEntropyH_primal() {
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

void nnDataGenerator::PrintTrainingData() {
    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        log->info( "{}, {}, {}", _uSol[idx_set][0], _alpha[idx_set][0], _hEntropy[idx_set] );
        logCSV->info( "{}, {}, {}", _uSol[idx_set][0], _alpha[idx_set][0], _hEntropy[idx_set] );
    }
}
