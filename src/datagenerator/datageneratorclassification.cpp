#include "datagenerator/datageneratorclassification.h"
#include "common/config.h"
#include "entropies/entropybase.h"
#include "optimizers/newtonoptimizer.h"
#include "quadratures/quadraturebase.h"
#include "spdlog/spdlog.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/sphericalbase.h"

#include <iostream>

DataGeneratorClassification::DataGeneratorClassification( Config* settings ) : DataGeneratorBase( settings ) {
    // Only 1D case right now
    _leftBound  = -5.0;
    _rightBound = 5.0;

    // Scale the quadrature weights
    _quadrature->ScaleWeights( _leftBound, _rightBound );
    _weights = _quadrature->GetWeights();
    _optimizer->ScaleQuadWeights( _leftBound, _rightBound );

    ComputeMoments();
    _uSol              = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _alpha             = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _pdfClassification = Vector( _setSize, 0.0 );

    double rho = 1.0;    // equal to the mass of the maxwellian
    double u   = 0.0;    // equal to the mean of the maxwellian
    double T   = 0.1;    // equal to the variance of the maxwellian

    _maxwellian = ComputeMaxwellian( rho, u, T );
}

DataGeneratorClassification::~DataGeneratorClassification() {}

void DataGeneratorClassification::ComputeTrainingData() {
    PrintLoadScreen();
    auto log = spdlog::get( "event" );
    log->info( "| Sample Lagrange multipliers." );
    SampleMultiplierAlpha();
    log->info( "| Multipliers sampled." );
    log->info( "| Reconstruct kinetic density functions." );
    log->info( "| Compute KL Divergence to Maxwellian with threshold XXX." );
    ClassifyDensity();
    log->info( "| Compute realizable problems." );
    ComputeRealizableSolution();
    log->info( "| Print Solution." );

    // --- Print everything ----
    PrintTrainingData();
}

void DataGeneratorClassification::ComputeMoments() {
    double my, phi;

    if( _settings->GetDim() == 1 ) {
        phi = 0;    // placeholder. will not be used

        for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
            my = _quadPointsSphere[idx_quad][0];

            // Scale my to match the trunctated velocity space
            my = ( _leftBound + _rightBound ) / 2.0 + my * ( _rightBound - _leftBound ) / 2.0;

            _momentBasis[idx_quad] = _basisGenerator->ComputeSphericalBasis( my, phi );

            // Correct the second moment with factor 0.5
            if( _nTotalEntries >= 3 ) {
                _momentBasis[idx_quad][2] = _momentBasis[idx_quad][2] * 0.5;
            }
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

    std::string uSolString = "";
    for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
        uSolString += "u_" + std::to_string( idx_sys ) + ",";
    }
    logCSV->info( uSolString + " equilibrium" );

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {

        std::stringstream streamU, streamEquilibrium;

        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            streamU << std::fixed << std::setprecision( 12 ) << _uSol[idx_set][idx_sys] << ",";
        }
        streamEquilibrium << std::fixed << std::setprecision( 1 ) << _pdfClassification[idx_set];

        std::string uSolString        = streamU.str();
        std::string equilibriumString = streamEquilibrium.str();

        logCSV->info( uSolString + equilibriumString );
    }
}
