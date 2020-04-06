#include "qldfesa.h"

QLDFESA::QLDFESA( unsigned order ) : QLookupQuadrature( order ){

    _availableOrders = {1, 2, 3}; //Available orders in lookuptable
    _nqByOrder = {32, 128, 512};
    _dataFiles = "../ext/sphericalquadpy/sphericalquadpy/ldfesa/data/";
    _dataFileSuffix = "_ldfesa.txt";

    SetName( ComputeName() );
    CheckOrder();         //Check if order is available
    SetNq( ComputeNq() ); //Set number of quadrature points

    SetPoints(ComputePoints());
    SetWeights(ComputeWeights());
    SetConnectivity( ComputeConnectivity() );
}

std::string QLDFESA::ComputeName() { return "LDFESA quadrature"; }

VectorVectorU QLDFESA::ComputeConnectivity() { //TODO
    VectorVectorU connectivity;
    return connectivity;
}
