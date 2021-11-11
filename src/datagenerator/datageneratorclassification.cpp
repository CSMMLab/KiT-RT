/*!
 * \file datageneratorclassification.h
 * \brief Class to generate data for the continuum breakdown classifier
 * \author S. Schotthoefer
 */

#include "datagenerator/datageneratorclassification.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "optimizers/newtonoptimizer.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/sphericalbase.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

#include "spdlog/spdlog.h"
#include <iostream>

DataGeneratorClassification::DataGeneratorClassification( Config* settings ) : DataGeneratorBase( settings ) {
    // Only 1D case right now
    _leftBound  = _settings->GetMinimalSamplingVelocity();
    _rightBound = _settings->GetMaximalSamplingVelocity();

    // Scale the quadrature weights
    _quadrature->ScalePointsAndWeights( _leftBound, _rightBound );
    _optimizer->ScaleQuadWeights( _leftBound, _rightBound );
    _weights          = _quadrature->GetWeights();
    _quadPointsSphere = _quadrature->GetPointsSphere();

    ComputeMoments();
    _uSol           = VectorVector();    // Not needed
    _alpha          = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _kineticDensity = VectorVector( _setSize, Vector( _nq, 0.0 ) );
}

DataGeneratorClassification::~DataGeneratorClassification() {}

void DataGeneratorClassification::ComputeTrainingData() {
    PrintLoadScreen();
    auto log = spdlog::get( "event" );
    log->info( "| Sample Lagrange multipliers." );
    SampleMultiplierAlpha();
    log->info( "| Multipliers sampled." );
    log->info( "| Reconstruct kinetic density functions." );
    ReconstructKineticDensity();
    log->info( "| Compute KL Divergence to Maxwellian with threshold XXX." );
    // ClassifyDensity();
    log->info( "| Compute realizable problems." );
    // ComputeRealizableSolution();
    log->info( "| Print Solution." );

    // --- Print everything ----
    PrintTrainingData();
}

void DataGeneratorClassification::SampleMultiplierAlpha() {

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
            double T   = idx_temp * dTemp;
            auto logg  = spdlog::get( "event" );
            logg->info( "| Sample Lagrange multipliers with mean based on Maxwellians with\n| Density =" + std::to_string( rho ) +
                        "\n| Bulk velocity =" + std::to_string( u ) + "\n| Temperature =" + std::to_string( T ) + "\n|" );

            Vector meanAlpha = Vector( 3, 0.0 );
            meanAlpha[0]     = log( rho / pow( 2 * M_PI, 3.0 / 2.0 ) ) - u * u / ( 2.0 * T );
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
                unsigned idx_set = idx_temp * nTemp + idx_loc;
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
                            std::cout << "here\n";
                        }
                        if( _alpha[idx_set][idx_sys] < meanAlpha[idx_sys] - maxAlphaValue ) {
                            _alpha[idx_set][idx_sys] = meanAlpha[idx_sys] - maxAlphaValue;
                            std::cout << "here2 \n";
                        }
                    }
                    // Compute rejection criteria
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
            meanAlpha[0]     = log( rho / pow( 2 * M_PI, 3.0 / 2.0 ) ) - u * u / ( 2.0 * T );
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
            for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
                Vector alphaRed = Vector( _nTotalEntries - 1, 0.0 );    // local reduced alpha

                bool accepted = false;
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

void DataGeneratorClassification::ReconstructKineticDensity() {
    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            _kineticDensity[idx_set][idx_quad] = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_set], _momentBasis[idx_quad] ) );
        }
    }
}

void DataGeneratorClassification::ComputeMoments() {
    double my, phi;

    if( _settings->GetDim() == 1 ) {
        phi = 0;    // placeholder. will not be used

        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            my = _quadPointsSphere[idx_quad][0];    // quad points are already scaled

            _momentBasis[idx_quad] = _basisGenerator->ComputeSphericalBasis( my, phi );

            // Correct the second moment with factor 0.5
            // if( _nTotalEntries >= 3 ) {
            //    _momentBasis[idx_quad][2] = _momentBasis[idx_quad][2] * 0.5;
            //}
        }
    }
    else if( _settings->GetDim() == 2 || _settings->GetDim() == 3 ) {
        ErrorMessages::Error( "Spatial Dimension not supported.", CURRENT_FUNCTION );

        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            my                     = _quadPointsSphere[idx_quad][0];
            phi                    = _quadPointsSphere[idx_quad][1];
            _momentBasis[idx_quad] = _basisGenerator->ComputeSphericalBasis( my, phi );
        }
    }
    else {
        ErrorMessages::Error( "Spatial Dimension not supported.", CURRENT_FUNCTION );
    }
}

void DataGeneratorClassification::ClassifyDensity() {
    double treshold = 1.0;    // dummy
    Vector currPdf  = Vector( _nq, 0.0 );
    for( unsigned idx_sol = 0; idx_sol < _setSize; idx_sol++ ) {
        _uSol[idx_sol] = 0;
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            // Make entropyReconstruction a member vector, s.t. it does not have to be re-evaluated in ConstructFlux
            currPdf[idx_quad] = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_sol], _momentBasis[idx_quad] ) );
        }
        double klDiv = ComputeKLDivergence( currPdf, _maxwellian );
        if( klDiv < treshold ) _pdfClassification[idx_sol] = 1.0;
    }
}

Vector DataGeneratorClassification::ComputeMaxwellian( double rho, double u, double T ) {
    // Only in 1D right now

    auto logCSV = spdlog::get( "tabular" );

    Vector maxwellian = Vector( _nq, 0.0 );
    double prefactor  = rho / sqrt( 2 * M_PI * T );
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        maxwellian[idx_quad] = prefactor * exp( -1 * ( ( _momentBasis[idx_quad][1] - u ) * ( _momentBasis[idx_quad][1] - u ) ) / ( 2 * T ) );
        std::cout << maxwellian[idx_quad] << std::endl;
    }

    // Compute the Moment of the maxwellian.
    Vector maxwellianMoment = Vector( _nTotalEntries, 0.0 );
    double moment0          = 0.0;
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        // std::cout << _momentBasis[idx_quad] << "\n";
        maxwellianMoment += _momentBasis[idx_quad] * _weights[idx_quad] * maxwellian[idx_quad];    //  _weights[idx_quad]
        moment0 += ( _weights[idx_quad] * maxwellian[idx_quad] );
    }
    // Compute the Lagrange multiplier of the maxwellian
    Vector maxwellianAlpha = Vector( _nTotalEntries, 0.0 );
    _optimizer->Solve( maxwellianAlpha, maxwellianMoment, _momentBasis );

    std::cout << "Maxwellian Moment:\n";
    std::cout << maxwellianMoment << std::endl;
    // std::cout << moment0 << std::endl;
    // std::cout << "Maxwellian Alpha:\n";
    std::cout << maxwellianAlpha << std::endl;

    maxwellianAlpha[1] *= 2;
    // For debugging, reconstruct the moments from Maxwellian alpha
    double entropyReconstruction = 0.0;
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        // Make entropyReconstruction a member vector, s.t. it does not have to be re-evaluated in ConstructFlux
        entropyReconstruction = _entropy->EntropyPrimeDual( blaze::dot( maxwellianAlpha, _momentBasis[idx_quad] ) );
        maxwellianMoment += _momentBasis[idx_quad] * ( _weights[idx_quad] * entropyReconstruction );
    }
    std::cout << "Reconstructed Maxwellian Moment:\n";
    std::cout << maxwellianMoment << std::endl;

    return maxwellian;
}

double DataGeneratorClassification::ComputeKLDivergence( Vector& f1, Vector& f2 ) {
    double sum = 0;
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        sum += f1[idx_quad] * log( f1[idx_quad] / f2[idx_quad] );
    }
    // std::cout << "KL-Divergence:" << sum << std::endl;
    return sum;
}

void DataGeneratorClassification::PrintTrainingData() {

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
