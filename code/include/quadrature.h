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
    virtual ~Quadrature(){}

    /*! @brief Gives the name of the quadrature rule
     *  @returns string name : Name of quadrature rule     */
    virtual std::string ComputeName()           = 0;

    /*! @brief Computes the number of gridpoints of the quadrature rule
     *  @returns unsigned nq : number of gridpoints of the quadrature rule     */
    virtual unsigned ComputeNq()                = 0;

    /*! @brief Computes the a vector (length: nq) of (coordinates of) gridpoints used for the quadrature rule
     *  @returns VectorVector coordinates : A Vector of coordinates the gridpoints.     */
    virtual VectorVector ComputePoints()        = 0;

    /*! @brief Computes the a vector (length: nq) of weights for the gridpoints. The indices match the gridpoints VectorVector.
     *  @returns Vector weights : A Vector of weights of the gridpoints.     */
    virtual Vector ComputeWeights()             = 0;

    /*! @brief TODO: How is connectivity defined?.
     *  @returns VectorVectorU connectivity : TODO */
    virtual VectorVectorU ComputeConnectivity() = 0;

    // Aux functions
    void PrintWeights();            /*! @brief prints: Weight vector */
    void PrintPoints();             /*! @brief prints: Point vectorVector */
    void PrintPointsAndWeights();   /*! @brief prints: Point vectorVector with corresponding weight vector */

     /*! @brief sums up all entries of the weight vector.
      *  @returns sum of all weights */
    double SumUpWeights();

     /*! @brief computes the weighted sum of a given flux function f over all quadrature points.
      *  @param double(f)( double x0, double x1, double x2 ) : flux function that depends on a three spatial dimensions.
      *  @returns double result: result of the quadrature rule */
    double Integrate( double( f )( double x0, double x1, double x2 ) );

    // Quadrature Hub
    /*! @brief Creates a quadrature rule with a given name and a given order.
     *  @param: std::string name: Name of the quadrature rule
     *  @param: unsigned order: Order of the quadrature rule
     *  @returns Quadrature* quadrature: returns pointer to instance of the given derived quadrature class */
    static Quadrature* CreateQuadrature( std::string name, unsigned order );

    // Setter
    inline void SetName( std::string name ) { _name = name; }           /*! @brief sets: name of the quadrature */
    inline void SetOrder( unsigned order ) { _order = order; }          /*! @brief sets: order of the quadrature */
    inline void SetNq( unsigned nq ) { _nq = nq; }                      /*! @brief sets: number of gridpoints of the quadrature */
    inline void SetPoints( VectorVector points ) { _points = points; }  /*! @brief sets: coordinates of gridpoints of the quadrature */
    inline void SetWeights( Vector weights ) { _weights = weights; }    /*! @brief sets: weights of gridpoints of the quadrature */
    inline void SetConnectivity( VectorVectorU connectivity ) { _connectivity = connectivity; } /*! @brief sets: connectivity vector */

    // Getter
    inline std::string GetName() { return _name; }      /*! @returns std::string _name:  name of the quadrature */
    inline unsigned GetOrder() { return _order; }       /*! @returns unsigned _order:  order of the quadrature */
    inline unsigned GetNq() { return _nq; }             /*! @returns unsigned _nq:  number of gridpoints of the quadrature */
    inline VectorVector GetPoints() { return _points; } /*! @returns VectorVector _points:  coordinates of gridpoints of the quadrature */
    inline Vector GetWeights() { return _weights; }     /*! @returns Vector _weights:  weights of gridpoints of the quadrature */
    inline VectorVectorU GetConnectivity() { return _connectivity; } /*! @returns VectorVectorU _connectivity:  connectivity of gridpoints of the quadrature */

  protected:
    std::string _name;          /*! @brief name of the quadrature */
    unsigned _order;            /*! @brief order of the quadrature */
    unsigned _nq;               /*! @brief number of gridpoints of the quadrature */
    VectorVector _points;       /*! @brief gridpoints of the quadrature */
    Vector _weights;            /*! @brief weights of the gridpoints of the quadrature */
    VectorVectorU _connectivity;/*! @brief connectivity of the gripoints of the quadrature */
};

#endif    // QUADRATURE_H
