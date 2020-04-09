#include "qlevelsymmetric.h"

QLevelSymmetric::QLevelSymmetric( unsigned order ) : QLookupQuadrature( order ){

    SetAvailOrders();
    SetDataInfo();

    SetName();
    CheckOrder(); //Check if order is available
    SetNq(); //Set number of quadrature points
    SetPointsAndWeights();
    SetConnectivity();
}

void QLevelSymmetric::SetAvailOrders() {
    _availableOrders = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20}; //Available orders in lookuptable
    _nqByOrder = {8, 24, 48, 80, 120, 168, 224, 288, 360, 432};
}

void QLevelSymmetric::SetDataInfo() {
    _dataFiles = "../ext/sphericalquadpy/sphericalquadpy/levelsymmetric/data/";
    _dataFileSuffix = "_levelsym.txt";
}

void QLevelSymmetric::SetConnectivity() { //TODO
    VectorVectorU connectivity;
    _connectivity = connectivity;
}
