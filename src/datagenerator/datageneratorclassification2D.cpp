/*!
 * \file datageneratorclassification2D.h
 * \brief Class to generate data for the continuum breakdown classifier in 2 spatial dimension
 * \author S. Schotthoefer
 */

#include "datagenerator/datageneratorclassification2D.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "optimizers/newtonoptimizer.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "velocitybasis/sphericalbase.hpp"

#include "spdlog/spdlog.h"
#include <iostream>

DataGeneratorClassification2D::DataGeneratorClassification2D( Config* settings ) : DataGeneratorClassification( settings ) {
    // ErrorMessages::Error( "2D Classification sampler is a work in progress\n", CURRENT_FUNCTION );
    ComputeMoments();
}

DataGeneratorClassification2D::~DataGeneratorClassification2D() {}

void DataGeneratorClassification2D::ComputeMoments() {
    double my, phi, r;
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        my                     = _quadPointsSphere[idx_quad][0];
        phi                    = _quadPointsSphere[idx_quad][1];
        r                      = _quadPointsSphere[idx_quad][2];
        _momentBasis[idx_quad] = _basisGenerator->ComputeSphericalBasis( my, phi, r );
    }
}

void DataGeneratorClassification2D::SampleMultiplierAlpha() {

    // Create mean alpha vector corresonding to a Maxwellian with rho=1, u=0, T=<sampled Temp>
    unsigned nTemp = _settings->GetNSamplingTemperatures();
    double tempMin = _settings->GetMinimalSamplingTemperature();
    double tempMax = _settings->GetMaximalSamplingTemperature();
    double dTemp   = ( tempMax - tempMin ) / double( nTemp );

    unsigned localSetSize = unsigned( _setSize / nTemp );

    if( !_settings->GetNormalizedSampling() ) {    // Not normalized
        for( unsigned idx_temp = 0; idx_temp < nTemp; idx_temp++ ) {
            if( idx_temp == nTemp - 1 ) {
                localSetSize = _setSize - ( nTemp - 1 ) * localSetSize;
            }
            double rho = 1.0;
            double u   = 0.0;    // for both dimension (both 0, so it's ok)
            double T   = idx_temp * dTemp + tempMin;

            auto logg = spdlog::get( "event" );
            logg->info( "| Sample Lagrange multipliers with mean based on Maxwellians with\n| Density =" + std::to_string( rho ) +
                        "\n| Bulk velocity =" + std::to_string( u ) + "\n| Temperature =" + std::to_string( T ) + "\n|" );

            Vector meanAlpha = Vector( 6, 0.0 );
            meanAlpha[0]     = log( rho / ( 2 * M_PI * T ) ) - u * u / ( 2.0 * T );
            meanAlpha[1]     = 2.0 * u / ( 2.0 * T );    // v_x
            meanAlpha[2]     = 2.0 * u / ( 2.0 * T );    // v_y
            meanAlpha[3]     = -1.0 / ( 2.0 * T );       // same for both directions
            meanAlpha[4]     = 0.0;                      // no mixed terms used
            meanAlpha[5]     = -1.0 / ( 2.0 * T );       // same for both directions

            logg->info( "| Lagrange multiplier means are:\n| Alpha0 =" + std::to_string( meanAlpha[0] ) +
                        "\n| Alpha1,Alpha2 =" + std::to_string( meanAlpha[1] ) + "\n| Alpha3 =" + std::to_string( meanAlpha[3] ) + "\n|" );

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
                        if( idx_sys == 0 ) {                            // v^0
                            _alpha[idx_set][idx_sys] = meanAlpha[0];    // fixed value
                        }
                        else if( idx_sys == 1 ) {                       // v_x
                            _alpha[idx_set][idx_sys] = meanAlpha[1];    // fixed value
                        }
                        else if( idx_sys == 2 ) {                       // v_y
                            _alpha[idx_set][idx_sys] = meanAlpha[2];    // fixed value
                        }
                        else if( idx_sys == 3 ) {                       // v_y^2
                            _alpha[idx_set][idx_sys] = meanAlpha[3];    // fixed value
                        }
                        else if( idx_sys == 4 ) {                                             // v_y*v_x
                            _alpha[idx_set][idx_sys] = distributionAlphaRest( generator );    // mixed term as disturbation
                        }
                        else if( idx_sys == 5 ) {                       // v_x^2
                            _alpha[idx_set][idx_sys] = meanAlpha[5];    // fixed value (same as v_y^2
                        }
                        else {
                            _alpha[idx_set][idx_sys] = distributionAlphaRest( generator );
                        }
                        if( _alpha[idx_set][idx_sys] > meanAlpha[idx_sys] + maxAlphaValue ) {
                            _alpha[idx_set][idx_sys] = meanAlpha[idx_sys] + maxAlphaValue;
                            std::cout << "lower bound touched\n";
                        }
                        if( _alpha[idx_set][idx_sys] < meanAlpha[idx_sys] - maxAlphaValue ) {
                            _alpha[idx_set][idx_sys] = meanAlpha[idx_sys] - maxAlphaValue;
                            std::cout << "upper bound touched \n";
                        }
                    }
                    // Compute rejection criteria
                    accepted = ComputeEVRejection( idx_set );
                }
            }
        }
    }
    else {    // Normalized
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

            Vector meanAlpha = Vector( 6, 0.0 );
            meanAlpha[0]     = log( rho / ( 2 * M_PI ) ) - u * u / ( 2.0 * T );
            meanAlpha[1]     = 2.0 * u / ( 2.0 * T );    // v_x
            meanAlpha[2]     = 2.0 * u / ( 2.0 * T );    // v_y
            meanAlpha[3]     = -1.0 / ( 2.0 * T );       // same for both directions
            meanAlpha[4]     = 0.0;                      // no mixed terms used
            meanAlpha[5]     = -1.0 / ( 2.0 * T );       // same for both directions

            logg->info( "| Lagrange multiplier means are:\n| Alpha1, Alpha2 =" + std::to_string( meanAlpha[1] ) +
                        "\n| Alpha3 =" + std::to_string( meanAlpha[3] ) + "\n|" );

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
                Vector alphaRed  = Vector( _nTotalEntries - 1, 0.0 );    // local reduced alpha

                bool accepted = false;
                while( !accepted ) {
                    // Sample random multivariate uniformly distributed alpha between minAlpha and MaxAlpha.
                    for( unsigned idx_sys = 1; idx_sys < _nTotalEntries; idx_sys++ ) {
                        if( idx_sys == 1 ) {                         // v_x
                            alphaRed[idx_sys - 1] = meanAlpha[1];    // fixed value
                        }
                        else if( idx_sys == 2 ) {                    // v_y
                            alphaRed[idx_sys - 1] = meanAlpha[2];    // fixed value
                        }
                        else if( idx_sys == 3 ) {                    // v_y^2
                            alphaRed[idx_sys - 1] = meanAlpha[3];    // fixed value
                        }
                        else if( idx_sys == 4 ) {                                          // v_y*v_x
                            alphaRed[idx_sys - 1] = distributionAlphaRest( generator );    // mixed term as disturbation
                        }
                        else if( idx_sys == 5 ) {                    // v_x^2
                            alphaRed[idx_sys - 1] = meanAlpha[5];    // fixed value (same as v_y^2
                        }
                        else {
                            alphaRed[idx_sys - 1] = distributionAlphaRest( generator );
                        }
                        if( alphaRed[idx_sys - 1] > meanAlpha[idx_sys] + maxAlphaValue ) {
                            alphaRed[idx_sys - 1] = meanAlpha[idx_sys] + maxAlphaValue;
                            // std::cout << "lower bound touched\n";
                        }
                        if( alphaRed[idx_sys - 1] < meanAlpha[idx_sys] - maxAlphaValue ) {
                            alphaRed[idx_sys - 1] = meanAlpha[idx_sys] - maxAlphaValue;
                            // std::cout << "upper bound touched \n";
                        }
                    }
                    // Compute alpha_0 = log(<exp(alpha m )>) // for maxwell boltzmann! only
                    double integral = 0.0;
                    // std::cout << alphaRed << "\n";
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
                    accepted = ComputeEVRejection( idx_set );
                }
            }
        }
    }
}

void DataGeneratorClassification2D::PrintTrainingData() {

    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );
    log->info( "---------------------- Data Generation Successful ------------------------" );

    std::stringstream quadPtsStream1, quadPtsStream2, quadPtsStream3, quadWeightsStream;
    for( unsigned idx_quad = 0; idx_quad < _nq - 1; idx_quad++ ) {
        quadPtsStream1 << std::fixed << std::setprecision( 12 ) << _quadPoints[idx_quad][0] << ",";
        quadPtsStream2 << std::fixed << std::setprecision( 12 ) << _quadPoints[idx_quad][1] << ",";
        quadPtsStream3 << std::fixed << std::setprecision( 12 ) << _quadPoints[idx_quad][2] << ",";
    }
    quadPtsStream1 << std::fixed << std::setprecision( 12 ) << _quadPoints[_nq - 1][0];
    quadPtsStream2 << std::fixed << std::setprecision( 12 ) << _quadPoints[_nq - 1][1];
    quadPtsStream3 << std::fixed << std::setprecision( 12 ) << _quadPoints[_nq - 1][2];
    std::string quadPtsString = quadPtsStream1.str();
    logCSV->info( quadPtsString );
    quadPtsString = quadPtsStream2.str();
    logCSV->info( quadPtsString );
    quadPtsString = quadPtsStream3.str();
    logCSV->info( quadPtsString );

    for( unsigned idx_quad = 0; idx_quad < _nq - 1; idx_quad++ ) {
        quadWeightsStream << std::fixed << std::setprecision( 12 ) << _weights[idx_quad] << ",";
    }
    quadWeightsStream << std::fixed << std::setprecision( 12 ) << _weights[_nq - 1];
    std::string quadWeightsString = quadWeightsStream.str();
    logCSV->info( quadWeightsString );

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {

        // Test the moments of the kinetic density
        // Vector u    = Vector( _nTotalEntries, 0.0 );
        // double tmp  = 0.0;
        // double recT = 0.0;
        // for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        //    u += _momentBasis[idx_quad] * _weights[idx_quad] * _kineticDensity[idx_set][idx_quad];
        //    tmp += _weights[idx_quad] * _kineticDensity[idx_set][idx_quad];
        //    recT += 0.5 * ( _momentBasis[idx_quad][1] * _momentBasis[idx_quad][1] + _momentBasis[idx_quad][2] * _momentBasis[idx_quad][2] ) *
        //            _weights[idx_quad] * _kineticDensity[idx_set][idx_quad];
        //}
        // std::cout << "alpha_" << idx_set << " : " << _alpha[idx_set] << "\n";
        // std::cout << "u_" << idx_set << " : " << u << "\n";
        // std::cout << "rho" << idx_set << " : " << tmp << "\n";
        // std::cout << "T" << idx_set << " : " << recT << "\n";

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
