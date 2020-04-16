#include "../../include/quadratures/qldfesa.h"

QLDFESA::QLDFESA( unsigned order ) : QLookupQuadrature( order ){

    SetAvailOrders();
    SetDataInfo();

    SetName();
    CheckOrder(); //Check if order is available
    SetNq(); //Set number of quadrature points
    SetPointsAndWeights();
    SetConnectivity();
}

void QLDFESA::SetAvailOrders() {
    _availableOrders = {1, 2, 3};
    _nqByOrder = {32, 128, 512};
}

void QLDFESA::SetDataInfo() {
    _dataFiles = "../ext/sphericalquadpy/sphericalquadpy/ldfesa/data/";
    _dataFileSuffix = "_ldfesa.txt";
}

void QLDFESA::SetConnectivity() { //TODO
    VectorVectorU connectivity;
    _connectivity = connectivity;
}
