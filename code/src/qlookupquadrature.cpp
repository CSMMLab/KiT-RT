#include <fstream>
#include <sstream>
#include "qlookupquadrature.h"

QLookupQuadrature::QLookupQuadrature( unsigned order ) : Quadrature( order ){
}

unsigned QLookupQuadrature::ComputeNq() {
    //find iterator of current order
    std::vector<unsigned>::iterator it = std::find(_availableOrders.begin(), _availableOrders.end(), _order);
    // Get index of element from iterator
    unsigned index = std::distance(_availableOrders.begin(), it);
    // Get _nq
    return _nqByOrder[index];
}

VectorVector QLookupQuadrature::ComputePoints() {
    std::string filename =  _dataFiles + std::to_string(_order) + _dataFileSuffix;

    Vector weights(_nq);
    VectorVector points( _nq );
    std::string line;
    unsigned count;

    ifstream in(filename); //give file to filestream
    for(unsigned idx_point = 0; idx_point < _nq ; idx_point++){

        count = 0;
        points[idx_point].resize(3);

        in >> line; //Get line of CSV file
        stringstream ss(line); //give line to strinstream

        for (double double_in; ss >> double_in;) { //parse line
            if(count <3) points[idx_point][count] = double_in;
            else weights[idx_point] = double_in;

            count ++;

            if (ss.peek() == ',') ss.ignore();
        }
    }
    in.close();

    return points;
}

Vector QLookupQuadrature::ComputeWeights() { //Basicially copied from above. Combine both functions for more efficiency?
    std::string filename =  _dataFiles + std::to_string(_order) + _dataFileSuffix;

    Vector weights(_nq);
    VectorVector points( _nq );
    std::string line;
    unsigned count;

    ifstream in(filename); //give file to filestream
    for(unsigned idx_point = 0; idx_point < _nq ; idx_point++){

        count = 0;
        points[idx_point].resize(3);

        in >> line; //Get line of CSV file
        stringstream ss(line); //give line to strinstream

        for (double double_in; ss >> double_in;) { //parse line
            if(count <3) points[idx_point][count] = double_in;
            else weights[idx_point] = double_in;

            count ++;

            if (ss.peek() == ',') ss.ignore();
        }
    }
    in.close();

    return weights;
}

bool QLookupQuadrature::CheckOrder(){
   std::vector<unsigned>::iterator it = std::find(_availableOrders.begin(), _availableOrders.end(), _order);

   if (it == _availableOrders.end()){
        cout << "ERROR! Order "<< _order << " not available. (Replace this error message by a proper exeption handler!)" << std::endl; //TODO: throw proper error!
        exit(1);
        return false;
    }
    std::cout << _name << " with order " << _order << " chosen.\n";
    return true;
}
