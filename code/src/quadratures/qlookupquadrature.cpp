#include "quadratures/qlookupquadrature.h"
#include "settings/globalconstants.h"    // for PI_NUMBER
#include <fstream>
#include <sstream>

QLookupQuadrature::QLookupQuadrature( unsigned order ) : QuadratureBase( order ) {}

void QLookupQuadrature::printAvailOrders() const {
    auto log                    = spdlog::get( "event" );
    std::string availableOrders = "";
    for( unsigned i = 0; i < _availableOrders.size() - 1; i++ ) {
        availableOrders += std::to_string( _availableOrders[i] ) + "|";
    }
    availableOrders += std::to_string( _availableOrders[_availableOrders.size() - 1] ) + ")";
    log->info( "Available orders: ({0})", availableOrders );
}

bool QLookupQuadrature::CheckOrder() {
    std::vector<unsigned>::iterator it = std::find( _availableOrders.begin(), _availableOrders.end(), _order );
    if( it == _availableOrders.end() ) {
        printAvailOrders();
        ErrorMessages::Error( "ERROR! Order " + std::to_string( _order ) + " for " + GetName() + " not available. ", CURRENT_FUNCTION );
    }
    return true;
}

void QLookupQuadrature::SetNq() {
    // find iterator of current order
    std::vector<unsigned>::iterator it = std::find( _availableOrders.begin(), _availableOrders.end(), _order );
    // Get index of element from iterator
    unsigned index = std::distance( _availableOrders.begin(), it );
    // Get _nq
    _nq = _nqByOrder[index];
}

void QLookupQuadrature::SetPointsAndWeights() {

    _weights.resize( _nq );
    _points.resize( _nq );

    double sumWeights = 0;

    std::string lookupTable = GetLookupTable();
    std::string line;
    unsigned count;

    std::stringstream in( lookupTable );

    for( unsigned idx_point = 0; idx_point < _nq; idx_point++ ) {

        count = 0;
        _points[idx_point].resize( 3 );
        line.clear();
        in >> line;    // Get line of lookupTable string

        if( line.empty() ) {
            ErrorMessages::Error( "Length of lookup table vector does not fit requested point vector length!", CURRENT_FUNCTION );
        }

        std::stringstream ss( line );    // give line to stringstream

        for( double double_in; ss >> double_in; ) {    // parse line
            if( count < 3 )
                _points[idx_point][count] = double_in;
            else {
                sumWeights += double_in;
                _weights[idx_point] = double_in;
            }
            count++;
            if( ss.peek() == ',' ) ss.ignore();
        }
    }

    // Correct the scaling of the weights
    for( unsigned idx = 0; idx < _nq; idx++ ) {
        _weights[idx] /= sumWeights;
        _weights[idx] *= 4.0 * PI_NUMBER;
    }
}
