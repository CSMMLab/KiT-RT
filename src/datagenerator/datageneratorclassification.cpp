#include "datagenerator/datageneratorclassification.h"

DataGeneratorClassification::DataGeneratorClassification( Config* settings ) : DataGeneratorBase( settings ) {
    _kineticDensity = VectorVector( _setSize, Vector( _nq, 0.0 ) );
    double rho      = 1.0;
    Vector u        = Vector( _nTotalEntries, 0.0 );
    double T        = 1.0;
    _maxwellian     = ComputeMaxwellian( rho, u, T );
}

DataGeneratorClassification::~DataGeneratorClassification() {}

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
    return sum;
}

void DataGeneratorClassification::ComputeTrainingData() {}

void DataGeneratorClassification::PrintTrainingData() {
    /*
    auto log    = spdlog::get( "event" );
    auto logCSV = spdlog::get( "tabular" );
    log->info( "---------------------- Data Generation Successful ------------------------" );

    std::string uSolString  = "";
    std::string alphaString = "";
    for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
        uSolString += "u_" + std::to_string( idx_sys ) + ",";
        alphaString += "alpha_" + std::to_string( idx_sys ) + ",";
    }
    // log->info( uSolString + alphaString + "h" );
    logCSV->info( uSolString + alphaString + "h" );

    for( unsigned idx_set = 0; idx_set < _setSize; idx_set++ ) {

        std::stringstream streamU, streamAlpha, streamH;

        for( unsigned idx_sys = 0; idx_sys < _nTotalEntries; idx_sys++ ) {
            streamU << std::fixed << std::setprecision( 12 ) << _uSol[idx_set][idx_sys] << ",";
            streamAlpha << std::fixed << std::setprecision( 12 ) << _alpha[idx_set][idx_sys] << ",";
        }
        streamH << std::fixed << std::setprecision( 12 ) << _hEntropy[idx_set];

        std::string uSolString  = streamU.str();
        std::string alphaString = streamAlpha.str();
        std::string hString     = streamH.str();

        // log->info(  uSolString + alphaString + hString  );
        logCSV->info( uSolString + alphaString + hString );
    }
    */
}
