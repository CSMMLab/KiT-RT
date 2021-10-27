//
// Created by chinsp on 27/10/21.
//

#ifndef QICOSAHEDRONTRIANG_H
#define QICOSAHEDRONTRIANG_H


#include "quadraturebase.h"

class QIcosahedronII : public QuadratureBase
{
public:
    QIcosahedronII( Config* settings );
    QIcosahedronII( unsigned quadOrder );
    virtual ~QIcosahedronII() {}

private:
    void SetName() override { _name = "Icosahedron quadrature."; }
    void SetNq() override;        /*!@brief Set number of quadrature points */
    bool CheckOrder();
    void SetConnectivity() override;

    void SetupTriangulation();    /*!@brief Setup of the Triangulation (face-wise) */
    void SetPointsAndWeights() override;

    // Helper
    void SetupPlatonicSolid();    /*!@brief Setup of the initial platonic solid (until now only Icosahedron supported) */
    VectorVector Interpolate( Vector, Vector, unsigned);  /*!@brief Interpolation between 2 points on unit sphere (slerp) */
    double GetArea( Vector, Vector, Vector);          /*!@brief Calculates the area of a spherical triangle */

    // Variables
private:
    VectorVectorU _faces;
    VectorVector _vertices;
    VectorVector _ptsTriang;
    VectorVectorU _triangles;
}

#endif //QICOSAHEDRONTRIANG_H
