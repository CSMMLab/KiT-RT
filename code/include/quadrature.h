#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <blaze/Math.h>
#include <math.h>
#include <random>
#include <string>

class Quadrature
{
  public:
    Quadrature( int order ) {
        SetOrder( order );
        SetName( ComputeName() );
        SetNq( ComputeNq() );
        SetPoints( ComputePoints() );
        SetWeights( ComputeWeights() );
        SetConnectivity( ComputeConnectivity() );
    };
    virtual ~Quadrature() = 0;

    // Virtual methods that every quadrature has to implement
    virtual std::string ComputeName();
    virtual int ComputeNq();
    virtual blaze::DynamicVector<blaze::DynamicVector<double>> ComputePoints();
    virtual blaze::DynamicVector<double> ComputeWeights();
    virtual blaze::DynamicVector<blaze::DynamicVector<int>> ComputeConnectivity();

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
