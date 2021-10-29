//
// Created by chinsp on 29/10/21.
//

#ifndef QICOSAHEDRON_H
#define QICOSAHEDRON_H

#include "quadraturebase.h"

class QIcosahedron : public QuadratureBase
{
public:
    QIcosahedron( Config* settings );
    QIcosahedron( unsigned quadOrder );
    virtual ~QIcosahedron() {}

private:
    void SetName() override { _name = "Icosahedron quadrature."; }
    void SetNq() override;
    bool CheckOrder();
    void SetConnectivity() override;

    void SetupTriangulation();
    void SetPointsAndWeights() override;

    /*! @brief Nr of quadraturepoints in refined quadrature */
    void SetNqRefined();
    /*! @brief Calculation of corresponding points and weights in refined quadrature /n
          Defines _pointsRefined, _weightsRefined, _refineVector */
    void Refine();
    /*! @brief Set values concerning refinement to original, if no refinement is used */
    void ResetValues();
    /*! @brief Determination of neighbour triangles for each quadrature point in COARSE quadrature \n
          Defines VectorVector _neighbours */
    void SetNeighbourConnectivity();

    // Helper
    /*! @brief Setup of the initial platonic solid (until now only Icosahedron supported) */
    void SetupPlatonicSolid();
    /*! @brief Interpolation between 2 points on unit sphere (slerp)
        @param Vector-a,b Points on unit sphere @param unsigned n Nr of Interpolation points
        @return VectorVector of Interpolation points
      */
    VectorVector Interpolate( Vector a, Vector b, unsigned n );
    /*! @brief Calculates the area of a spherical triangle
         @param Vector-a,b,c Points on unit sphere
         @return double Area of spherical triangle
     */

    double GetArea( Vector a, Vector b, Vector c);

    // Variables
private:
    VectorVectorU _faces;
    VectorVector _vertices;
    VectorVector _ptsTriang;
    VectorVectorU _triangles;

};

#endif //QICOSAHEDRON_H
