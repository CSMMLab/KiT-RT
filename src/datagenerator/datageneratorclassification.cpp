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
    _quadrature->ScalePointsAndWeights( _rightBound );
    _optimizer->ScaleQuadWeights( _rightBound );
    _weights = _quadrature->GetWeights();
    // std::cout << "sum of weights: " << _quadrature->SumUpWeights() << "\n";
    _quadPointsSphere = _quadrature->GetPointsSphere();
    _quadPoints       = _quadrature->GetPoints();

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

void DataGeneratorClassification::ReconstructKineticDensity() {
    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {
        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            _kineticDensity[idx_set][idx_quad] = _entropy->EntropyPrimeDual( blaze::dot( _alpha[idx_set], _momentBasis[idx_quad] ) );
        }
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
