/*!
 * \file datageneratorbase.cpp
 * \brief Class to generate data for the neural entropy closure
 * \author S. Schotthoefer
 */

#include "datagenerator/datageneratorbase.hpp"
#include "common/config.hpp"
#include "datagenerator/datageneratorclassification1D.hpp"
#include "datagenerator/datageneratorclassification2D.hpp"
#include "datagenerator/datageneratorregression1D.hpp"
#include "datagenerator/datageneratorregression2D.hpp"
#include "datagenerator/datageneratorregression3D.hpp"
#include "entropies/entropybase.hpp"
#include "optimizers/newtonoptimizer.hpp"
#include "optimizers/partregularizednewtonoptimizer.hpp"
#include "optimizers/reducednewtonoptimizer.hpp"
#include "optimizers/reducedpartregularizednewtonoptimizer.hpp"
#include "optimizers/regularizednewtonoptimizer.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "velocitybasis/sphericalbase.hpp"

//#include <chrono>
#include "spdlog/spdlog.h"
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
    _reducedSampling = false;
    switch( settings->GetOptimizerName() ) {
        case NEWTON: _optimizer = new NewtonOptimizer( settings ); break;
        case REGULARIZED_NEWTON: _optimizer = new RegularizedNewtonOptimizer( settings ); break;
        case PART_REGULARIZED_NEWTON: _optimizer = new PartRegularizedNewtonOptimizer( settings ); break;
        case REDUCED_NEWTON:
            _optimizer       = new ReducedNewtonOptimizer( settings );
            _reducedSampling = true;
            break;
        case REDUCED_PART_REGULARIZED_NEWTON:
            _optimizer       = new ReducedPartRegularizedNewtonOptimizer( settings );
            _reducedSampling = true;
            break;
        default: ErrorMessages::Error( "Optimizer choice not feasible for datagenerator.", CURRENT_FUNCTION ); break;
    }
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
            case 1: return new DataGeneratorRegression1D( settings );
            case 2: return new DataGeneratorRegression2D( settings );
            case 3: return new DataGeneratorRegression3D( settings );
            default: ErrorMessages::Error( "Sampling for more than 3 dimensions is not yet supported.", CURRENT_FUNCTION );
        }
    }
    else if( settings->GetSamplerName() == CLASSIFICATION_SAMPLER ) {
        switch( settings->GetDim() ) {
            case 1: return new DataGeneratorClassification1D( settings );
            case 2: return new DataGeneratorClassification2D( settings );
            default: ErrorMessages::Error( "Sampling for more than 3 dimensions is not yet supported.", CURRENT_FUNCTION );
        }
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
        double mean   = 0.0;
        double stddev = maxAlphaValue / 3.0;
        std::normal_distribution<double> distribution_normal( mean, stddev );

        // Can be parallelized, but check if there is a race condition with datagenerator
        for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
            Vector alphaRed = Vector( _nTotalEntries - 1, 0.0 );    // local reduced alpha

            bool accepted = false;
            bool normAccepted = false;
            while( !accepted ) {
                // Sample random multivariate uniformly distributed alpha between minAlpha and MaxAlpha.
                for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                    if( _settings->GetUniformSamlping() )
                        alphaRed[idx_sys - 1] = distribution( generator );
                    else {
                        alphaRed[idx_sys - 1] = distribution_normal( generator );
                        if( alphaRed[idx_sys - 1] > maxAlphaValue ) alphaRed[idx_sys - 1] = maxAlphaValue;
                        if( alphaRed[idx_sys - 1] < -1 * maxAlphaValue ) alphaRed[idx_sys - 1] = -1 * maxAlphaValue;
                    }
                }
                normAccepted = true;
                if( _settings->GetUniformSamlping() ) {
                    if( norm( alphaRed ) > maxAlphaValue ) normAccepted = false;
                }
                // Compute alpha_0 = log(<exp(alpha m )>) // for maxwell boltzmann! only
                double integral = 0.0;
                // Integrate <eta'_*(alpha*m)>
                for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
                    integral += _entropy->EntropyPrimeDual( dot( alphaRed, momentsRed[idx_quad] ) ) * _weights[idx_quad];
                }
                _alpha[idx_set][0] = -log( integral );    // log trafo

                // Assemble complete alpha (normalized)
                for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                    _alpha[idx_set][idx_sys] = alphaRed[idx_sys - 1];
                }

                // Compute rejection criteria
                accepted = normAccepted;
                if( normAccepted ) {
                    if( _reducedSampling ) {
                        accepted = ComputeReducedEVRejection( momentsRed, alphaRed );
                    }
                    else {
                        accepted = ComputeEVRejection( idx_set );
                    }
                }
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
    double minAlphaValue = -1 * maxAlphaValue;
    // double maxAlphaValue = maxAlphaValue;
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
        double _gridSize = floor( pow( (double)_setSize, 1 / (double)( _nTotalEntries - 1 ) ) );
        double dalpha    = ( maxAlphaValue - minAlphaValue ) / (double)_gridSize;
        unsigned count   = 0;

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
        return false;
    }
    return true;
}

bool DataGeneratorBase::ComputeReducedEVRejection( VectorVector& redMomentBasis, Vector& redAlpha ) {
    Matrix hessian = Matrix( _nTotalEntries - 1, _nTotalEntries - 1, 0.0 );
    _optimizer->ComputeHessian( redAlpha, redMomentBasis, hessian );
    SymMatrix hessianSym( hessian );    // Bad solution, rewrite with less memory need
    Vector ew = Vector( _nTotalEntries - 1, 0.0 );
    eigen( hessianSym, ew );
    if( min( ew ) < _settings->GetMinimalEVBound() ) {
        return false;
    }
    return true;
}

void DataGeneratorBase::ComputeRealizableSolution() {

    if( _reducedSampling ) {
        VectorVector momentsRed = VectorVector( _nq, Vector( _nTotalEntries - 1, 0.0 ) );

        for( unsigned idx_nq = 0; idx_nq < _nq; idx_nq++ ) {    // copy (reduced) moments
            for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                momentsRed[idx_nq][idx_sys - 1] = _momentBasis[idx_nq][idx_sys];
            }
        }
#pragma omp parallel for schedule( guided )
        for( unsigned idx_sol = 0; idx_sol < _setSize; idx_sol++ ) {
            Vector uSolRed( _nTotalEntries - 1, 0.0 );
            Vector alphaRed( _nTotalEntries - 1, 0.0 );
            for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                alphaRed[idx_sys - 1] = _alpha[idx_sol][idx_sys];
            }
            // std::cout << alphaRed << std::endl;
            // std::cout << _alpha[idx_sol] << std::endl;
            _optimizer->ReconstructMoments( uSolRed, alphaRed, momentsRed );
            // std::cout << uSolRed << std::endl;

            for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                _uSol[idx_sol][idx_sys] = uSolRed[idx_sys - 1];
            }
            _uSol[idx_sol][0] = 1.0;
        }
    }
    else {
#pragma omp parallel for schedule( guided )
        for( unsigned idx_sol = 0; idx_sol < _setSize; idx_sol++ ) {
            _optimizer->ReconstructMoments( _uSol[idx_sol], _alpha[idx_sol], _momentBasis );
        }
    }

    // TextProcessingToolbox::PrintVectorVector( _uSol );
}
