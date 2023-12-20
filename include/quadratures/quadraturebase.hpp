/*! @file: quadraturebase.h
 *  @brief Base class for all quadrature rules in KiT-RT
 *  @author: S. Schotthöfer
 */

#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "common/globalconstants.hpp"
#include "common/typedef.hpp"

class Config;

class QuadratureBase
{
  public:
    /*! @brief Constructor using settings class. This is the recommended constructor.
     *  @param settings Settings class storing all important options.
     */
    QuadratureBase( Config* settings );
    /*! @brief Constructor using directly the order of the quadrature. Not applicable for GaussLegendre, that need additional options.
                It sets member _settings = nulltpr.*/
    QuadratureBase( unsigned order );
    virtual ~QuadratureBase() {}

    // Aux functions
    void PrintWeights();          /*!< @brief prints: Weight vector */
    void PrintPoints();           /*!< @brief prints: Point vectorVector */
    void PrintPointsAndWeights(); /*!< @brief prints: Point vectorVector with corresponding weight vector */

    /*! @brief sums up all entries of the weight vector.
     *  @returns sum of all weights */
    double SumUpWeights();

    /*! @brief Integrates f(x,y,z) with the quadrature.
     *  @param f density function that depends on a three spatial dimensions.
     *  @returns result of the quadrature rule */
    virtual double Integrate( double ( *f )( double, double, double ) );

    /*! @brief Integrates f(x,y,z) with the quadrature.
     *  @param f density function that depends on a spherical coordinates.
     *  @returns result of the quadrature rule */
    virtual double IntegrateSpherical( double ( *f )( double, double ) );

    /*! @brief Integrates vector valued f(x,y,z) with the quadrature. Each dimension is integrated by itself.
     *  @param f density function that depends on a three spatial dimensions.
     *  @param len lenght of vector
     *  @returns result of the quadrature rule (vector valued) */
    virtual std::vector<double> Integrate( std::vector<double> ( *f )( double, double, double ), unsigned len );

    // Quadrature Hub
    /*! @brief Creates a quadrature rule with a given name and a given order.
     *  @param  settings  Settings to handle quadrature options
     *  @returns Quadrature* quadrature: returns pointer to instance of the given derived quadrature class */
    static QuadratureBase* Create( Config* settings );

    /*! @brief Creates a quadrature rule with a given name and a given order.
     *  @param name name of quadrature as enum
     *  @param quadOrder order of quadrature
     *  @returns pointer to instance of the given derived quadrature class */
    static QuadratureBase* Create( QUAD_NAME name, unsigned quadOrder );

    // Getter
    inline std::string GetName() const { return _name; }                  /*! @returns std::string _name:  name of the quadrature */
    inline unsigned GetOrder() const { return _order; }                   /*! @returns unsigned _order:  order of the quadrature */
    inline unsigned GetNq() const { return _nq; }                         /*! @returns unsigned _nq:  number of gridpoints of the quadrature */
    inline const VectorVector& GetPoints() const { return _pointsKarth; } /*! @returns VectorVector _points:  coordinates of the quad pts */
    inline const VectorVector& GetPointsSphere() const {
        return _pointsSphere;
    } /*! @returns VectorVector _pointsSphere:  "---- " in spherical coordinates (my, phi)*/
    inline const Vector& GetWeights() const { return _weights; } /*! @returns Vector _weights:  weights of gridpoints of the quadrature */

    /*! @returns VectorVectorU _connectivity:  connectivity of gridpoints of the quadrature */
    inline const VectorVectorU& GetConnectivity() const { return _connectivity; }
    /*!< @brief Returns approved Dimensions for this quadrature */
    inline const std::vector<unsigned short> GetSupportedDims() const { return _supportedDimensions; }
    /*!< @brief Scales the quadrature weights according to the intervall [-velocityScaling , velocityScaling] in 1D
                Scales the radius of the velocity sphere by velocityScaling in 2D and 3D */
    virtual void ScalePointsAndWeights( double velocityScaling );

  protected:
    // Setter
    inline void SetOrder( unsigned order ) { _order = order; } /*! @brief sets: order of the quadrature */
    virtual void SetName()         = 0;                        /*!< @brief Sets: name of the quadrature */
    virtual void SetNq()           = 0;                        /*!< @brief sets: number of gridpoints of the quadrature */
    virtual void SetConnectivity() = 0;                        /*!< @brief sets: Connectivity Adjacency Matrix as VektorVektor*/

    /*! @brief Computes the a vector (length: nq) of (coordinates of) gridpoints used for the quadrature rule.
     *         Computes the a vector (length: nq) of weights for the gridpoints. The indices match the gridpoints VectorVector.
     *         Sets computed values for _points and _weights. */
    virtual void SetPointsAndWeights() = 0;

    // Member variables
    Config* _settings;                                /*!< @brief pointer to settings class that manages the solver */
    std::string _name;                                /*!< @brief name of the quadrature */
    unsigned _order;                                  /*!< @brief order of the quadrature */
    unsigned _nq;                                     /*!< @brief number of gridpoints of the quadrature */
    VectorVector _pointsKarth;                        /*!< @brief gridpoints of the quadrature */
    VectorVector _pointsSphere;                       /*!< @brief (my,phi,r)gridpoints of the quadrature in spherical cordinates */
    Vector _weights;                                  /*!< @brief weights of the gridpoints of the quadrature */
    VectorVectorU _connectivity;                      /*!< @brief connectivity of the gripoints of the quadrature */
    std::vector<unsigned short> _supportedDimensions; /*!< @brief number of spatial dimensions, for which the quadrature is build */
};

#endif    // QUADRATURE_H
