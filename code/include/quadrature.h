#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <blaze/Math.h>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
class Quadrature
{
  public:
    Quadrature( int order );
    virtual ~Quadrature(){};

    virtual std::string ComputeName()                                             = 0;
    virtual int ComputeNq()                                                       = 0;
    virtual blaze::DynamicVector<blaze::DynamicVector<double>> ComputePoints()    = 0;
    virtual blaze::DynamicVector<double> ComputeWeights()                         = 0;
    virtual blaze::DynamicVector<blaze::DynamicVector<int>> ComputeConnectivity() = 0;

    // Aux functions
    void PrintWeights();
    void PrintPoints();
    double SumUpWeights();
    

    // Quadrature Hub
    static Quadrature* CreateQuadrature( std::string name, int order );

    // Setter
    void SetName( std::string name ) { _name = name; };
    void SetOrder( int order ) { _order = order; };
    void SetNq( int nq ) { _nq = nq; };
    void SetPoints( blaze::DynamicVector<blaze::DynamicVector<double>> points ) { _points = points; };
    void SetWeights( blaze::DynamicVector<double> weights ) { _weights = weights; };
    void SetConnectivity( blaze::DynamicVector<blaze::DynamicVector<int>> connectivity ) { _connectivity = connectivity; };

    // Getter
    std::string GetName() { return _name; };
    int GetOrder() { return _order; };
    int GetNq() { return _nq; };
    blaze::DynamicVector<blaze::DynamicVector<double>> GetPoints() { return _points; };
    blaze::DynamicVector<double> GetWeights() { return _weights; };
    blaze::DynamicVector<blaze::DynamicVector<int>> GetConnectivity() { return _connectivity; };

  protected:
    std::string _name;
    int _order;
    int _nq;
    blaze::DynamicVector<blaze::DynamicVector<double>> _points;
    blaze::DynamicVector<double> _weights;
    blaze::DynamicVector<blaze::DynamicVector<int>> _connectivity;
};

#endif    // QUADRATURE_H
