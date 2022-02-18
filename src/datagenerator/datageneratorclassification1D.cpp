/*!
 * \file datageneratorclassification1D.h
 * \brief Class to generate data for the continuum breakdown classifier in 1 spatial dimension
 * \author S. Schotthoefer
 */

#include "datagenerator/datageneratorclassification1D.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "optimizers/newtonoptimizer.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "velocitybasis/sphericalbase.hpp"

#include "spdlog/spdlog.h"
#include <iostream>

DataGeneratorClassification1D::DataGeneratorClassification1D( Config* settings ) : DataGeneratorClassification( settings ) { ComputeMoments(); }

DataGeneratorClassification1D::~DataGeneratorClassification1D() {}

void DataGeneratorClassification1D::ComputeMoments() {
    double my, phi;
    phi = 0;    // placeholder. will not be used
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my                     = _quadPointsSphere[idx_quad][0];    // quad points are already scaled
        _momentBasis[idx_quad] = _basisGenerator->ComputeSphericalBasis( my, phi );
    }
}

void DataGeneratorClassification1D::SampleMultiplierAlpha() {

    // Create mean alpha vector corresonding to a Maxwellian with rho=1, u=0, T=<sampled Temp>
    unsigned nTemp = _settings->GetNSamplingTemperatures();
    double tempMin = _settings->GetMinimalSamplingTemperature();
    double tempMax = _settings->GetMaximalSamplingTemperature();
    double dTemp   = ( tempMax - tempMin ) / double( nTemp );

    unsigned localSetSize = unsigned( _setSize / nTemp );

    if( !_settings->GetNormalizedSampling() ) {
        for( unsigned idx_temp = 0; idx_temp < nTemp; idx_temp++ ) {
            if( idx_temp == nTemp - 1 ) {
                localSetSize = _setSize - ( nTemp - 1 ) * localSetSize;
            }
            double rho = 1.0;
            double u   = 0.0;
            double T   = idx_temp * dTemp + tempMin;
            auto logg  = spdlog::get( "event" );
            logg->info( "| Sample Lagrange multipliers with mean based on Maxwellians with\n| Density =" + std::to_string( rho ) +
                        "\n| Bulk velocity =" + std::to_string( u ) + "\n| Temperature =" + std::to_string( T ) + "\n|" );

            Vector meanAlpha = Vector( 3, 0.0 );
            meanAlpha[0]     = log( pow( rho / ( 2 * M_PI * T ), 1.0 / 2.0 ) ) - u * u / ( 2.0 * T );
            meanAlpha[1]     = 2.0 * u / ( 2.0 * T );
            meanAlpha[2]     = -1.0 / ( 2.0 * T );

            logg->info( "| Lagrange multiplier means are:\n| Alpha0 =" + std::to_string( meanAlpha[0] ) +
                        "\n| Alpha1 =" + std::to_string( meanAlpha[1] ) + "\n| Alpha2 =" + std::to_string( meanAlpha[2] ) + "\n|" );

            double maxAlphaValue = _settings->GetAlphaSamplingBound();
            double stddev        = maxAlphaValue / 3.0;

            std::default_random_engine generator;
            std::normal_distribution<double> distributionAlpha0( meanAlpha[0], stddev );
            std::normal_distribution<double> distributionAlpha1( meanAlpha[1], stddev );
            std::normal_distribution<double> distributionAlpha2( meanAlpha[2], stddev );
            std::normal_distribution<double> distributionAlphaRest( 0.0, stddev );

            // Can be parallelized, but check if there is a race condition with datagenerator
            for( unsigned idx_loc = 0; idx_loc < localSetSize; idx_loc++ ) {
                unsigned idx_set = idx_temp * localSetSize + idx_loc;
                bool accepted    = false;
                while( !accepted ) {
                    // Sample random multivariate uniformly distributed alpha between minAlpha and MaxAlpha.
                    for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
                        if( idx_sys == 0 ) {
                            _alpha[idx_set][idx_sys] = meanAlpha[0];    // distributionAlpha0( generator );
                        }
                        else if( idx_sys == 1 ) {
                            _alpha[idx_set][idx_sys] = meanAlpha[1];    // distributionAlpha1( generator );
                        }
                        else if( idx_sys == 2 ) {
                            _alpha[idx_set][idx_sys] = meanAlpha[2];    // distributionAlpha2( generator );
                        }
                        else {
                            _alpha[idx_set][idx_sys] = distributionAlphaRest( generator );
                        }
                        if( _alpha[idx_set][idx_sys] > meanAlpha[idx_sys] + maxAlphaValue ) {
                            _alpha[idx_set][idx_sys] = meanAlpha[idx_sys] + maxAlphaValue;
                            std::cout << "upper bound touched\n";
                        }
                        if( _alpha[idx_set][idx_sys] < meanAlpha[idx_sys] - maxAlphaValue ) {
                            _alpha[idx_set][idx_sys] = meanAlpha[idx_sys] - maxAlphaValue;
                            std::cout << "lower bound touched \n";
                        }
                    }
                    //  Compute rejection criteria
                    accepted = ComputeEVRejection( idx_set );
                }
            }
        }
    }
    else {
        for( unsigned idx_temp = 0; idx_temp < nTemp; idx_temp++ ) {
            if( idx_temp == nTemp - 1 ) {
                localSetSize = _setSize - ( nTemp - 1 ) * localSetSize;
            }
            double rho = 1.0;
            double u   = 0.0;
            double T   = idx_temp * dTemp + tempMin;
            auto logg  = spdlog::get( "event" );
            logg->info( "| Sample Lagrange multipliers with mean based on Maxwellians with\n| Density =" + std::to_string( rho ) +
                        "\n| Bulk velocity =" + std::to_string( u ) + "\n| Temperature =" + std::to_string( T ) + "\n|" );

            Vector meanAlpha = Vector( 3, 0.0 );
            meanAlpha[0]     = log( rho / pow( 2 * M_PI * T, 1.0 / 2.0 ) ) - u * u / ( 2.0 * T );
            meanAlpha[1]     = 2.0 * u / ( 2.0 * T );
            meanAlpha[2]     = -1.0 / ( 2.0 * T );

            logg->info( "| Lagrange multiplier means are:\n| Alpha1 =" + std::to_string( meanAlpha[1] ) +
                        "\n| Alpha2 =" + std::to_string( meanAlpha[2] ) + "\n|" );

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
            double maxAlphaValue = _settings->GetAlphaSamplingBound();
            double stddev        = maxAlphaValue / 3.0;
            std::normal_distribution<double> distributionAlpha1( meanAlpha[1], stddev );
            std::normal_distribution<double> distributionAlpha2( meanAlpha[2], stddev );
            std::normal_distribution<double> distributionAlphaRest( 0.0, stddev );

            // Can be parallelized, but check if there is a race condition with datagenerator
            for( unsigned idx_loc = 0; idx_loc < localSetSize; idx_loc++ ) {
                unsigned idx_set = idx_temp * localSetSize + idx_loc;

                Vector alphaRed = Vector( _nTotalEntries - 1, 0.0 );    // local reduced alpha
                bool accepted   = false;

                while( !accepted ) {
                    // Sample random multivariate uniformly distributed alpha between minAlpha and MaxAlpha.
                    for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {

                        if( idx_sys == 1 ) {
                            alphaRed[idx_sys - 1] = meanAlpha[1];    // distributionAlpha1( generator );
                        }
                        else if( idx_sys == 2 ) {
                            alphaRed[idx_sys - 1] = meanAlpha[2];    // distributionAlpha2( generator );
                        }
                        else {
                            alphaRed[idx_sys - 1] = distributionAlphaRest( generator );
                        }
                    }
                    // Compute alpha_0 = log(<exp(alpha m )>) // for maxwell boltzmann! only
                    double integral = 0.0;

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
    }
}

void DataGeneratorClassification1D::PrintTrainingData() {

    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );
    log->info( "---------------------- Data Generation Successful ------------------------" );

    std::stringstream quadPtsStream, quadWeightsStream;
    for( unsigned idx_quad = 0; idx_quad < _nq - 1; idx_quad++ ) {
        quadPtsStream << std::fixed << std::setprecision( 12 ) << _quadPointsSphere[idx_quad][0] << ",";
    }
    quadPtsStream << std::fixed << std::setprecision( 12 ) << _quadPointsSphere[_nq - 1][0];
    std::string quadPtsString = quadPtsStream.str();
    logCSV->info( quadPtsString );

    for( unsigned idx_quad = 0; idx_quad < _nq - 1; idx_quad++ ) {
        quadWeightsStream << std::fixed << std::setprecision( 12 ) << _weights[idx_quad] << ",";
    }
    quadWeightsStream << std::fixed << std::setprecision( 12 ) << _weights[_nq - 1];
    std::string quadWeightsString = quadWeightsStream.str();
    logCSV->info( quadWeightsString );

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {

        std::stringstream streamDensity;
        for( unsigned idx_quad = 0; idx_quad < _nq - 1; idx_quad++ ) {
            streamDensity << std::fixed << std::setprecision( 12 ) << _kineticDensity[idx_set][idx_quad] << ",";
        }
        streamDensity << std::fixed << std::setprecision( 12 ) << _kineticDensity[idx_set][_nq - 1];

        std::string densityString = streamDensity.str();

        logCSV->info( densityString );
    }
    log->info( "------------------------- Data printed to file ---------------------------" );
}
