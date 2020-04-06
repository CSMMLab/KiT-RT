#include "qlevelsymmetric.h"

QLevelSymmetric::QLevelSymmetric( unsigned order ) : QLookupQuadrature( order ){

    _availableOrders = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20}; //Available orders in lookuptable
    _nqByOrder = {8, 24, 48, 80, 120, 168, 224, 288, 360, 432};
    _dataFiles = "../ext/sphericalquadpy/sphericalquadpy/levelsymmetric/data/";
    _dataFileSuffix = "_levelsym.txt";

    SetName( ComputeName() );
    CheckOrder(); //Check if order is available
    SetNq( ComputeNq() ); //Set number of quadrature points

    SetPoints(ComputePoints());
    SetWeights(ComputeWeights());
    SetConnectivity( ComputeConnectivity() );
}

std::string QLevelSymmetric::ComputeName() { return "Level Symmetric quadrature"; }

VectorVectorU QLevelSymmetric::ComputeConnectivity() { //TODO
    VectorVectorU connectivity;
    return connectivity;
}
