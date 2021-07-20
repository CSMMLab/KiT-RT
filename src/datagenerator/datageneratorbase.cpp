/*!
 * \file datageneratorbase.cpp
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#include "datagenerator/datageneratorbase.h"
#include "common/config.h"
#include "datagenerator/datagenerator1D.h"
#include "datagenerator/datagenerator2D.h"
#include "datagenerator/datagenerator3D.h"
#include "datagenerator/datageneratorclassification.h"
#include "entropies/entropybase.h"
#include "optimizers/newtonoptimizer.h"
#include "quadratures/quadraturebase.h"
#include "spdlog/spdlog.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/sphericalbase.h"
#include "toolboxes/textprocessingtoolbox.h"

#include <chrono>
#include <iomanip>
#include <math.h>
#include <omp.h>
#include <sstream>

DataGeneratorBase::DataGeneratorBase( Config* settings ) {
    _settings = settings;
    _setSize  = settings->GetTrainingDataSetSize();

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
    _basisGenerator = SphericalBase::Create( _settings );

    _nTotalEntries = _basisGenerator->GetBasisSize();

    _momentBasis = VectorVector( _nq, Vector( _nTotalEntries, 0.0 ) );

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

    if( settings->GetSamplerName() == REGRESSION_SAMPLER ) {
        switch( settings->GetDim() ) {
            case 1: return new DataGenerator1D( settings );
            case 2: return new DataGenerator2D( settings );
            case 3: return new DataGenerator3D( settings );
            default: ErrorMessages::Error( "Sampling for more than 3 dimensions is not yet supported.", CURRENT_FUNCTION );
        }
    }
    else if( settings->GetSamplerName() == CLASSIFICATION_SAMPLER ) {
        return new DataGeneratorClassification( settings );
    }
    return nullptr;
}

void DataGeneratorBase::SampleMultiplierAlpha() {
    double maxAlphaValue = _settings->GetAlphaSamplingBound();
    // Rejection Sampling based on smallest EV of H

    if( _settings->GetNormalizedSampling() ) {
        if( _maxPolyDegree == 0 ) {
            ErrorMessages::Error( "Normalized sampling not meaningful for M0 closure", CURRENT_FUNCTION );
        }

        VectorVector momentsRed = VectorVector( _nq, Vector( _nTotalEntries - 1, 0.0 ) );

        for( unsigned idx_nq = 0; idx_nq < _nq; idx_nq++ ) {    // copy (reduced) moments
            for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                momentsRed[idx_nq][idx_sys - 1] = _momentBasis[idx_nq][idx_sys];
            }
        }

        // Create generator
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution( -1 * maxAlphaValue, maxAlphaValue );

        // Can be parallelized, but check if there is a race condition with datagenerator
        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            Vector alphaRed = Vector( _nTotalEntries - 1, 0.0 );    // local reduced alpha

            bool accepted = false;
            while( !accepted ) {
                // Sample random multivariate uniformly distributed alpha between minAlpha and MaxAlpha.
                for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                    alphaRed[idx_sys - 1] = distribution( generator );
                }
                // Compute alpha_0 = log(<exp(alpha m )>) // for maxwell boltzmann! only
                double integral = 0.0;
                // Integrate (eta(eta'_*(alpha*m))
                for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                    integral += _entropy->EntropyPrimeDual( dot( alphaRed, momentsRed[idx_quad] ) ) * _weights[idx_quad];
                }
                _alpha[idx_set][0] = -log( integral );    // log trafo

                // Assemble complete alpha (normalized)
                for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                    _alpha[idx_set][idx_sys] = alphaRed[idx_sys - 1];
                }

                // Compute rejection criteria
                accepted = ComputeEVRejection( idx_set );
            }
        }
    }
    else {
        // non normalized sampling
        // TODO
        ErrorMessages::Error( "Non-Normalized Alpha Sampling is not yet implemented.", CURRENT_FUNCTION );
    }

    // Cubeoid
    /*
    double minAlphaValue = -20;
    double maxAlphaValue = 20;
    if( _settings->GetNormalizedSampling() ) {
        // compute reduced version of alpha and m
        if( _maxPolyDegree == 0 ) {
            ErrorMessages::Error( "Normalized sampling not meaningful for M0 closure", CURRENT_FUNCTION );
        }

        VectorVector alphaRed   = VectorVector( _setSize, Vector( _nTotalEntries - 1, 0.0 ) );
        VectorVector momentsRed = VectorVector( _nq, Vector( _nTotalEntries - 1, 0.0 ) );

        for( unsigned idx_nq = 0; idx_nq < _nq; idx_nq++ ) {    // copy (reduced) moments
            for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                momentsRed[idx_nq][idx_sys - 1] = _momentBasis[idx_nq][idx_sys];
            }
        }

        // Sample alphaRed as uniform grid from [minAlphaValue, maxAlphaValue], then compute alpha_0 s.t. u_0 = 1
        double dalpha  = ( maxAlphaValue - minAlphaValue ) / (double)_gridSize;
        unsigned count = 0;

        switch( _nTotalEntries - 1 ) {
            case 1:
                for( unsigned idx_set = 0; idx_set < _gridSize; idx_set++ ) {
                    alphaRed[idx_set][0] = minAlphaValue + idx_set * dalpha;
                }
                break;
            case 2:
                count = 0;
                for( unsigned i1 = 0; i1 < _gridSize; i1++ ) {
                    double alpha0 = minAlphaValue + i1 * dalpha;
                    for( unsigned i2 = 0; i2 < _gridSize; i2++ ) {
                        alphaRed[count][0] = alpha0;
                        alphaRed[count][1] = minAlphaValue + i2 * dalpha;
                        count++;
                    }
                }
                break;
            case 3:
                count = 0;
                for( unsigned i1 = 0; i1 < _gridSize; i1++ ) {
                    double alpha0 = minAlphaValue + i1 * dalpha;
                    for( unsigned i2 = 0; i2 < _gridSize; i2++ ) {
                        double alpha1 = minAlphaValue + i2 * dalpha;
                        for( unsigned i3 = 0; i3 < _gridSize; i3++ ) {
                            alphaRed[count][0] = alpha0;
                            alphaRed[count][1] = alpha1;
                            alphaRed[count][2] = minAlphaValue + i3 * dalpha;
                            count++;
                        }
                    }
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
    */
}

void DataGeneratorBase::PrintLoadScreen() {
    auto log = spdlog::get( "event" );
    log->info( "------------------------ Data Generation Starts --------------------------" );
    log->info( "| Generating {} datapoints.", _setSize );
}

bool DataGeneratorBase::ComputeEVRejection( unsigned idx_set ) {

    Matrix hessian = Matrix( _nTotalEntries, _nTotalEntries, 0.0 );
    _optimizer->ComputeHessian( _alpha[idx_set], _momentBasis, hessian );

    SymMatrix hessianSym( hessian );    // Bad solution, rewrite with less memory need
    Vector ew = Vector( _nTotalEntries, 0.0 );
    eigen( hessianSym, ew );
    if( min( ew ) < _settings->GetMinimalEVBound() ) {
        std::cout << "Sampling not accepted with EV:" << min( ew ) << std::endl;
        // std::cout << "Current minimal accepted EV:" << _settings->GetMinimalEVBound() << std::endl;
        return false;
    }
    return true;
}
