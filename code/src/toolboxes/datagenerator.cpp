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
#include <omp.h>

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

nnDataGenerator::~nnDataGenerator() {
    delete _quadrature;
    delete _entropy;
}

void nnDataGenerator::computeTrainingData() {
    // Prototype: Only for _LMaxDegree == 1
    // Prototype: u is sampled from [0,100]

    // --- sample u ---
    SampleSolutionU();

    PrintLoadScreen();

    // ---- Check realizability ---
    CheckRealizability();

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
    double du0 = _settings->GetMaxValFirstMoment() / (double)_setSize;

    // different processes for different
    if( _LMaxDegree == 0 ) {
        // --- sample u in order 0 ---
        // u_0 = <1*psi>

        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            _uSol[idx_set][0] = du0 * idx_set;
        }
    }
    if( _LMaxDegree == 1 ) {
        // Sample points on unit sphere.
        QuadratureBase* quad = QuadratureBase::Create( _settings );
        VectorVector qpoints = quad->GetPoints();    // carthesian coordinates.
        unsigned long nq     = (unsigned long)quad->GetNq();

        // Allocate memory.
        unsigned long setSizeU0 = _setSize;
        unsigned long setSizeU1 = nq * setSizeU0 * ( setSizeU0 - 1 ) / 2;
        _setSize                = setSizeU1;

        // REFACTOR THIS
        _uSol     = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
        _alpha    = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
        _hEntropy = std::vector<double>( _setSize, 0.0 );

        // --- sample u in order 0 ---
        // u_0 = <1*psi>

#pragma omp parallel for schedule( guided )
        for( unsigned long idx_set = 0; idx_set < setSizeU0; idx_set++ ) {

            unsigned long outerIdx = ( idx_set - 1 ) * ( idx_set ) / 2;    // sum over all radii up to current
            outerIdx *= nq;                                                // in each radius step, use all quad points

            //  --- sample u in order 1 ---
            /* order 1 has 3 elements. (omega_x, omega_y,  omega_z) = omega, let u_1 = (u_x, u_y, u_z) = <omega*psi>
             * Condition u_0 >= norm(u_1)   */
            double radiusU0 = du0 * ( idx_set + 1 ) * ( 1 + 2 * _settings->GetBoundaryDistanceRealizableSet() );

            unsigned long localIdx = 0;
            double radius          = 0.0;
            unsigned long innerIdx = 0;
            // loop over all radii
            for( unsigned long idx_subset = 0; idx_subset < idx_set; idx_subset++ ) {
                radius   = du0 * ( idx_subset + 1 );    // dont use radius 0 ==> shift by one
                localIdx = outerIdx + idx_subset * nq;

                for( unsigned long quad_idx = 0; quad_idx < nq; quad_idx++ ) {
                    innerIdx = localIdx + quad_idx;    // gives the global index

                    _uSol[innerIdx][0] = radiusU0;    // Prevent 1st order moments to be too close to the boundary
                    // scale quadpoints with radius
                    _uSol[innerIdx][1] = radius * qpoints[quad_idx][0];
                    _uSol[innerIdx][2] = radius * qpoints[quad_idx][1];
                    _uSol[innerIdx][3] = radius * qpoints[quad_idx][2];
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
    log->info( "---------------------- Data Generation Successful ------------------------" );

    std::string uSolString  = "";
    std::string alphaString = "";
    for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
        uSolString += "u_" + std::to_string( idx_sys ) + " , ";
        alphaString += "alpha_" + std::to_string( idx_sys ) + " , ";
    }
    // log->info( uSolString + alphaString + " h" );
    logCSV->info( uSolString + alphaString + " h" );

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        std::string uSolString  = "";
        std::string alphaString = "";
        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            uSolString += std::to_string( _uSol[idx_set][idx_sys] ) + " , ";
            alphaString += std::to_string( _alpha[idx_set][idx_sys] ) + " , ";
        }
        // log->info( uSolString + alphaString + " {}", _hEntropy[idx_set] );
        logCSV->info( uSolString + alphaString + " {}", _hEntropy[idx_set] );
    }
}

void nnDataGenerator::CheckRealizability() {
    double epsilon = _settings->GetBoundaryDistanceRealizableSet();
    if( _LMaxDegree == 1 ) {
        double normU1 = 0.0;
        Vector u1( 3, 0.0 );
        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            if( _uSol[idx_set][0] < epsilon ) {
                if( std::abs( _uSol[idx_set][1] ) > 0 || std::abs( _uSol[idx_set][2] ) > 0 || std::abs( _uSol[idx_set][3] ) > 0 ) {
                    ErrorMessages::Error( "Moment not realizable [code 0].", CURRENT_FUNCTION );
                }
            }
            else {
                u1     = { _uSol[idx_set][1], _uSol[idx_set][2], _uSol[idx_set][3] };
                normU1 = norm( u1 );
                if( normU1 / _uSol[idx_set][0] > 1 - 0.99 * epsilon ) {
                    // std::cout << "normU1 / _uSol[" << idx_set << "][0]: " << normU1 / _uSol[idx_set][0] << "\n";
                    // std::cout << "normU1: " << normU1 << " | _uSol[idx_set][0] " << _uSol[idx_set][0] << "\n";
                    ErrorMessages::Error( "Moment to close to boundary of realizable set [code 1].\nBoundary ratio: " +
                                              std::to_string( normU1 / _uSol[idx_set][0] ),
                                          CURRENT_FUNCTION );
                }
                if( normU1 / _uSol[idx_set][0] <= 0 /*+ 0.5 * epsilon*/ ) {
                    // std::cout << "_uSol" << _uSol[idx_set][1] << " | " << _uSol[idx_set][2] << " | " << _uSol[idx_set][3] << " \n";
                    // std::cout << "normU1 / _uSol[" << idx_set << "][0]: " << normU1 / _uSol[idx_set][0] << "\n";
                    // std::cout << "normU1: " << normU1 << " | _uSol[idx_set][0] " << _uSol[idx_set][0] << "\n";
                    ErrorMessages::Error( "Moment to close to boundary of realizable set [code 2].\nBoundary ratio: " +
                                              std::to_string( normU1 / _uSol[idx_set][0] ),
                                          CURRENT_FUNCTION );
                }
            }
        }
    }
}

void nnDataGenerator::PrintLoadScreen() {
    auto log = spdlog::get( "event" );
    log->info( "------------------------ Data Generation Starts --------------------------" );
    log->info( "| Generating {} datapoints.", _setSize );
}
