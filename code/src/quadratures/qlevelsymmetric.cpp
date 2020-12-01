#include "quadratures/qlevelsymmetric.h"
//#include "quadratures/lookuptable_levelsymmetric.h"
#include "toolboxes/errormessages.h"

QLevelSymmetric::QLevelSymmetric( Config* settings ) : QLookupQuadrature( settings ) {
    SetAvailOrders();
    SetName();
    CheckOrder();    // Check if order is available
    SetNq();         // Set number of quadrature points
    SetPointsAndWeights();
    SetConnectivity();
}

QLevelSymmetric::QLevelSymmetric( unsigned quadOrder ) : QLookupQuadrature( quadOrder ) {
    SetAvailOrders();
    SetName();
    CheckOrder();    // Check if order is available
    SetNq();         // Set number of quadrature points
    SetPointsAndWeights();
    SetConnectivity();
}

void QLevelSymmetric::SetAvailOrders() {
    _availableOrders = { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 };    // Available orders in lookuptable
    _nqByOrder       = { 8, 24, 48, 80, 120, 168, 224, 288, 360, 432 };
}

void QLevelSymmetric::SetConnectivity() {    // TODO
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

std::string QLevelSymmetric::GetLookupTable() {

    unsigned char* lookupTable = nullptr;

    switch( _order ) {
        // case 2: lookupTable = __2_levelsym_txt; break;
        // case 4: lookupTable = __4_levelsym_txt; break;
        // case 6: lookupTable = __6_levelsym_txt; break;
        // case 8: lookupTable = __8_levelsym_txt; break;
        // case 10: lookupTable = __10_levelsym_txt; break;
        // case 12: lookupTable = __12_levelsym_txt; break;
        // case 14: lookupTable = __14_levelsym_txt; break;
        // case 16: lookupTable = __16_levelsym_txt; break;
        // case 18: lookupTable = __18_levelsym_txt; break;
        // case 20: lookupTable = __20_levelsym_txt; break;
        default: ErrorMessages::Error( "Invalid quadrature order chosen!", CURRENT_FUNCTION );
    }

    std::string lookupTableString( reinterpret_cast<char*>( lookupTable ) );

    return lookupTableString;
}
