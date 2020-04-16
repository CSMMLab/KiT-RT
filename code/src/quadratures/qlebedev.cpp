#include "../../include/quadratures/qlebedev.h"

QLebedev::QLebedev( unsigned order ) : QLookupQuadrature( order ){

    SetAvailOrders();
    SetDataInfo();

    SetName();
    CheckOrder(); //Check if order is available
    SetNq(); //Set number of quadrature points
    SetPointsAndWeights();
    SetConnectivity();
}

void QLebedev::SetAvailOrders()  {
    _availableOrders = {3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53,
                        59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131}; //Available orders in lookuptable
    _nqByOrder = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590,
                  770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};
}

void QLebedev::SetDataInfo() {
    _dataFiles = "../ext/sphericalquadpy/sphericalquadpy/lebedev/data/";
   _dataFileSuffix = "_lebedev.txt";
}

void QLebedev::SetConnectivity() { //TODO
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

