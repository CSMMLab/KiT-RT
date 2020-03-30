#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "typedef.h"
#include <iostream>
#include <string>
using namespace std;
class Quadrature
{
  public:
    Quadrature( unsigned order );
    virtual ~Quadrature(){};

    virtual std::string ComputeName()           = 0;
    virtual int ComputeNq()                     = 0;
    virtual VectorVector ComputePoints()        = 0;
    virtual Vector ComputeWeights()             = 0;
    virtual VectorVectorU ComputeConnectivity() = 0;

    // Aux functions
    void PrintWeights();
    void PrintPoints();
    void PrintPointsAndWeights();
    double SumUpWeights();

    // Quadrature Hub
    static Quadrature* CreateQuadrature( std::string name, unsigned order );

    // Setter
    void SetName( std::string name ) { _name = name; };
    void SetOrder( unsigned order ) { _order = order; };
    void SetNq( unsigned nq ) { _nq = nq; };
    void SetPoints( VectorVector points ) { _points = points; };
    void SetWeights( Vector weights ) { _weights = weights; };
    void SetConnectivity( VectorVectorU connectivity ) { _connectivity = connectivity; };

    // Getter
    std::string GetName() { return _name; };
    unsigned GetOrder() { return _order; };
    unsigned GetNq() { return _nq; };
    VectorVector GetPoints() { return _points; };
    Vector GetWeights() { return _weights; };
    VectorVectorU GetConnectivity() { return _connectivity; };

  protected:
    std::string _name;
    unsigned _order;
    unsigned _nq;
    VectorVector _points;
    Vector _weights;
    VectorVectorU _connectivity;
};

#endif    // QUADRATURE_H
