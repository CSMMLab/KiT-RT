#include "quadrature.h"
#include <iostream>

int main( int argc, char** argv ) {
    std::cout << "Hello world!" << std::endl;
    Quadrature* Q = GetQuadrature( "montecarlo", 10 );
    Q->PrintWeights();
    return EXIT_SUCCESS;
}
