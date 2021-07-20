#include "datagenerator/datageneratorclassification.h"
#include "common/config.h"
#include "entropies/entropybase.h"
#include "spdlog/spdlog.h"
#include "toolboxes/errormessages.h"
#include "toolboxes/sphericalbase.h"

#include <iostream>

DataGeneratorClassification::DataGeneratorClassification( Config* settings ) : DataGeneratorBase( settings ) {
    ComputeMoments();
    _uSol              = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _alpha             = VectorVector( _setSize, Vector( _nTotalEntries, 0.0 ) );
    _pdfClassification = Vector( _setSize, 0.0 );

    double rho = 1.0;
    Vector u   = Vector( _nTotalEntries, 0.0 );
    double T   = 1.0;

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
            my                     = _quadPointsSphere[idx_quad][0];
            _momentBasis[idx_quad] = _basisGenerator->ComputeSphericalBasis( my, phi );
        }
    }
    else if( _settings->GetDim() == 2 || _settings->GetDim() == 3 ) {
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

Vector DataGeneratorClassification::ComputeMaxwellian( double rho, Vector u, double T ) {
    Vector maxwellian = Vector( _nq, 0.0 );
    double prefactor  = rho / pow( 2 * M_PI * T, 2.0 / 3.0 );
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        maxwellian[idx_quad] = prefactor * exp( -1 * dot( _momentBasis[idx_quad] - u, _momentBasis[idx_quad] - u ) / ( 2 * T ) );
    }
    return maxwellian;
}

double DataGeneratorClassification::ComputeKLDivergence( Vector& f1, Vector& f2 ) {
    double sum = 0;
    for( unsigned idx_quad = 0; idx_quad < _nq; idx_quad++ ) {
        sum += f1[idx_quad] * log( f1[idx_quad] / f2[idx_quad] );
    }
    std::cout << "KL-Divergence:" << sum << std::endl;
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
