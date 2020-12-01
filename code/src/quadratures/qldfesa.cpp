#include "quadratures/qldfesa.h"
//#include "quadratures/lookuptable_ldfesa.h"
#include "toolboxes/errormessages.h"

QLDFESA::QLDFESA( Config* settings ) : QLookupQuadrature( settings ) {
    SetAvailOrders();
    SetName();
    CheckOrder();    // Check if order is available
    SetNq();         // Set number of quadrature points
    SetPointsAndWeights();
    SetConnectivity();
}
QLDFESA::QLDFESA( unsigned quadOrder ) : QLookupQuadrature( quadOrder ) {
    SetAvailOrders();
    SetName();
    CheckOrder();    // Check if order is available
    SetNq();         // Set number of quadrature points
    SetPointsAndWeights();
    SetConnectivity();
}

void QLDFESA::SetAvailOrders() {
    _availableOrders = { 1, 2, 3 };
    _nqByOrder       = { 32, 128, 512 };
}

void QLDFESA::SetConnectivity() {    // TODO
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

std::string QLDFESA::GetLookupTable() {

    unsigned char* lookupTable = nullptr;

    switch( _order ) {
        // case 1: lookupTable = __1_ldfesa_txt; break;
        // case 2: lookupTable = __2_ldfesa_txt; break;
        // case 3: lookupTable = __3_ldfesa_txt; break;
        default: ErrorMessages::Error( "Invalid quadrature order chosen!", CURRENT_FUNCTION );
    }

    std::string lookupTableString( reinterpret_cast<char*>( lookupTable ) );

    return lookupTableString;
}
