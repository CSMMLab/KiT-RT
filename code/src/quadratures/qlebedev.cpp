#include "quadratures/qlebedev.h"
#include "quadratures/lookuptable_lebedev.h"
#include "toolboxes/errormessages.h"

QLebedev::QLebedev( Config* settings ) : QLookupQuadrature( settings ) {

    SetAvailOrders();

    SetName();
    CheckOrder();    // Check if order is available
    SetNq();         // Set number of quadrature points
    SetPointsAndWeights();
    SetConnectivity();
}

void QLebedev::SetAvailOrders() {
    _availableOrders = { 3,  5,  7,  9,  11, 13, 15, 17, 19, 21, 23,  25,  27,  29,  31,  35,
                         41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131 };    // Available orders in lookuptable
    _nqByOrder       = { 6,   14,  26,  38,   50,   74,   86,   110,  146,  170,  194,  230,  266,  302,  350,  434,
                   590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 };
}

void QLebedev::SetConnectivity() {    // TODO
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

std::string QLebedev::GetLookupTable() {

    unsigned char* lookupTable = nullptr;

    switch( _order ) {
        case 3: lookupTable = __3_lebedev_txt; break;
        case 5: lookupTable = __5_lebedev_txt; break;
        case 7: lookupTable = __7_lebedev_txt; break;
        case 9: lookupTable = __9_lebedev_txt; break;
        case 11: lookupTable = __11_lebedev_txt; break;
        case 13: lookupTable = __13_lebedev_txt; break;
        case 15: lookupTable = __15_lebedev_txt; break;
        case 17: lookupTable = __17_lebedev_txt; break;
        case 19: lookupTable = __19_lebedev_txt; break;
        case 21: lookupTable = __21_lebedev_txt; break;
        case 23: lookupTable = __23_lebedev_txt; break;
        case 25: lookupTable = __25_lebedev_txt; break;
        case 27: lookupTable = __27_lebedev_txt; break;
        case 29: lookupTable = __29_lebedev_txt; break;
        case 31: lookupTable = __31_lebedev_txt; break;
        case 35: lookupTable = __35_lebedev_txt; break;
        case 41: lookupTable = __41_lebedev_txt; break;
        case 47: lookupTable = __47_lebedev_txt; break;
        case 53: lookupTable = __53_lebedev_txt; break;
        case 59: lookupTable = __59_lebedev_txt; break;
        case 65: lookupTable = __65_lebedev_txt; break;
        case 71: lookupTable = __71_lebedev_txt; break;
        case 77: lookupTable = __77_lebedev_txt; break;
        case 83: lookupTable = __83_lebedev_txt; break;
        case 89: lookupTable = __89_lebedev_txt; break;
        case 95: lookupTable = __95_lebedev_txt; break;
        case 101: lookupTable = __101_lebedev_txt; break;
        case 107: lookupTable = __107_lebedev_txt; break;
        case 113: lookupTable = __113_lebedev_txt; break;
        case 119: lookupTable = __119_lebedev_txt; break;
        case 125: lookupTable = __125_lebedev_txt; break;
        case 131: lookupTable = __131_lebedev_txt; break;
        default: ErrorMessages::Error( "Invalid quadrature order chosen!", CURRENT_FUNCTION );
    }

    std::string lookupTableString( reinterpret_cast<char*>( lookupTable ) );

    return lookupTableString;
}
