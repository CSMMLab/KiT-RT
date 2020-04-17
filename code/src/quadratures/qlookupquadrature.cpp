#include <fstream>
#include <sstream>
#include "../../include/quadratures/qlookupquadrature.h"
#include "option_structure.h" // for PI_NUMBER

QLookupQuadrature::QLookupQuadrature( unsigned order ) : Quadrature( order ){
}

void QLookupQuadrature::printAvailOrders() const{
    std::cout << "Available orders: (";
    for(unsigned i = 0; i < _availableOrders.size()-1 ; i++){
        std::cout << _availableOrders[i] << "|";
    }
    std::cout << _availableOrders[_availableOrders.size()-1] << ")\n";
}

bool QLookupQuadrature::CheckOrder(){
   std::vector<unsigned>::iterator it = std::find(_availableOrders.begin(), _availableOrders.end(), _order);

   if (it == _availableOrders.end()){
        std::cerr << "ERROR! Order "<< _order << " for " << GetName() << " not available. " << std::endl;
        printAvailOrders();
        exit(EXIT_FAILURE);
    }
    return true;
}

void QLookupQuadrature::SetNq() {
    //find iterator of current order
    std::vector<unsigned>::iterator it = std::find(_availableOrders.begin(), _availableOrders.end(), _order);
    // Get index of element from iterator
    unsigned index = std::distance(_availableOrders.begin(), it);
    // Get _nq
    _nq = _nqByOrder[index];
}

void QLookupQuadrature::SetPointsAndWeights() {

    _weights.resize(_nq);
    _points.resize(_nq);

    double sumWeights = 0;

    std::string filename =  _dataFiles + std::to_string(_order) + _dataFileSuffix;
    std::string line;
    unsigned count;

    std::ifstream in(filename); //give file to filestream

    // check, if file is empty
    if (in.peek() == std::ifstream::traits_type::eof()){
        std::cerr << "ERROR! Lookup table not found. Please check if externals are loaded correctly. \n" << std::endl;
        exit(EXIT_FAILURE);
    }

    for(unsigned idx_point = 0; idx_point < _nq ; idx_point++){

        count = 0;
        _points[idx_point].resize(3);

        in >> line; //Get line of CSV file
        std::stringstream ss(line); //give line to stringstream

        for (double double_in; ss >> double_in;) { //parse line
            if(count <3) _points[idx_point][count] = double_in;
            else{
                sumWeights += double_in;
                _weights[idx_point] = double_in;
            }
            count ++;
            if (ss.peek() == ',') ss.ignore();
        }
    }
    in.close();

    //Correct the scaling of the weights
    for (unsigned idx =0; idx < _nq; idx ++){
        _weights[idx] /= sumWeights;
        _weights[idx] *= 4.0 * PI_NUMBER;
    }
}
