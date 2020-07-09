#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "settings/globalconstants.h"
#include "settings/typedef.h"
#include "toolboxes/errormessages.h"

class QuadratureBase
{
  public:
    QuadratureBase( unsigned order );
    virtual ~QuadratureBase() {}

    // Aux functions
    void PrintWeights();          /*! @brief prints: Weight vector */
    void PrintPoints();           /*! @brief prints: Point vectorVector */
    void PrintPointsAndWeights(); /*! @brief prints: Point vectorVector with corresponding weight vector */

    /*! @brief sums up all entries of the weight vector.
     *  @returns sum of all weights */
    double SumUpWeights();

    /*! @brief Integrates f(x,y,z) with the quadrature.
     *  @param double(f)( double x0, double x1, double x2 ) : density function that depends on a three spatial dimensions.
     *  @returns double result: result of the quadrature rule */
    double Integrate( double( f )( double x0, double x1, double x2 ) );

    /*! @brief Integrates vector valued f(x,y,z) with the quadrature. Each dimension is integrated by itself.
     *  @param : double(f)( double x0, double x1, double x2 ) : density function that depends on a three spatial dimensions.
     *  @param :  len : lenght of vector
     *  @returns double result: result of the quadrature rule (vector valued) */
    std::vector<double> Integrate( std::vector<double>( f )( double x0, double x1, double x2 ), unsigned len );

    // Quadrature Hub
    /*! @brief Creates a quadrature rule with a given name and a given order.
     *  @param: std::string name: Name of the quadrature rule
     *  @param: unsigned order: Order of the quadrature rule
     *  @returns Quadrature* quadrature: returns pointer to instance of the given derived quadrature class */
    static QuadratureBase* CreateQuadrature( QUAD_NAME name, unsigned order );

    // Getter
    inline std::string GetName() const { return _name; }      /*! @returns std::string _name:  name of the quadrature */
    inline unsigned GetOrder() const { return _order; }       /*! @returns unsigned _order:  order of the quadrature */
    inline unsigned GetNq() const { return _nq; }             /*! @returns unsigned _nq:  number of gridpoints of the quadrature */
    inline VectorVector GetPoints() const { return _points; } /*! @returns VectorVector _points:  coordinates of gridpoints of the quadrature */
    virtual VectorVector GetPointsSphere() const;             /*! @returns VectorVector _pointsSphere:  "---- " in spherical coordinates (my, phi)*/
    inline Vector GetWeights() const { return _weights; }     /*! @returns Vector _weights:  weights of gridpoints of the quadrature */
    inline VectorVectorU GetConnectivity() const {
        return _connectivity;
    } /*! @returns VectorVectorU _connectivity:  connectivity of gridpoints of the quadrature */

  protected:
    // Setter
    inline void SetOrder( unsigned order ) { _order = order; } /*! @brief sets: order of the quadrature */
    virtual void SetName()         = 0;                        /*! @brief Sets: name of the quadrature */
    virtual void SetNq()           = 0;                        /*! @brief sets: number of gridpoints of the quadrature */
    virtual void SetConnectivity() = 0;                        /*! @brief sets: Connectivity Adjacency Matrix as VektorVektor*/

    /*! @brief Computes the a vector (length: nq) of (coordinates of) gridpoints used for the quadrature rule.
     *         Computes the a vector (length: nq) of weights for the gridpoints. The indices match the gridpoints VectorVector.
     *         Sets computed values for _points and _weights. */
    virtual void SetPointsAndWeights() = 0;

    // Member variables
    // TODO Config* _settings;           /*! @brief pointer to settings class that manages the solver */
    std::string _name;           /*! @brief name of the quadrature */
    unsigned _order;             /*! @brief order of the quadrature */
    unsigned _nq;                /*! @brief number of gridpoints of the quadrature */
    VectorVector _points;        /*! @brief gridpoints of the quadrature */
    Vector _weights;             /*! @brief weights of the gridpoints of the quadrature */
    VectorVectorU _connectivity; /*! @brief connectivity of the gripoints of the quadrature */

    VectorVector _pointsSphere; /*! @brief gridpoints of the quadrature in spherical cordinates */
};

#endif    // QUADRATURE_H
