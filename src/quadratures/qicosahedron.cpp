//
// Created by chinsp on 29/10/21.
//

#include "../../include/quadratures/qicosahedron.h"
#include "../../include/common/config.h"
#include "../../include/toolboxes/errormessages.h"

#include <stdio.h>
#include <iostream>
#include <time.h>

/*
Icosahedron Quadrature
Triangulation of each face:
    ---decompose each triangle from lower order in 4 new triangles---
    Order 1 -> 1 triangle per face
    Order 2 -> 4 triangles per face
    Order 3 -> 16 triangles per face ...
Uses midpoints of subtriangles as Quadrature Points
Uses areas of subtriangles as Quadrature Weights
Refinement ist defined
*/


QIcosahedron::QIcosahedron( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    CheckOrder();
    SetNq();
    SetupTriangulation();
    SetPointsAndWeights();
    SetConnectivity();
    // Refinement
    SetNqRefined();
    Refine();
    SetNeighbourConnectivity();
    std::cout<< "Quad ready"<<std::endl;
}

QIcosahedron::QIcosahedron( unsigned order ) : QuadratureBase( order ) {
    SetName();
    CheckOrder();
    SetNq();
    SetupTriangulation();
    SetPointsAndWeights();
    SetConnectivity();
}

void QIcosahedron::SetNq() {
    _nq = 20 * pow( 4, _order - 1  );
}

void QIcosahedron::SetNqRefined(){
    if ( _settings->GetQuadOrderFine() == _order + 2 ){
        _nqRefined = 20 * pow( 4, _order + 1);
    }
    else if ( _settings->GetQuadOrderFine() == _order +1 ){
        _nqRefined = 20 * pow( 4, _order);
    }
    else if (_settings->GetQuadOrderFine() == _order ){
        _nqRefined = _nq;
    }
    else{
        ErrorMessages::Error("Please define QUAD_ORDER_FINE as QUAD_ORDER + 1 or QUAD_ORDER + 2", "Refinement of Icosahedron Quadrature.");
    }
}

void QIcosahedron::SetupTriangulation() {

    unsigned orderTriang = pow( 2, _order-1) + 1;
    unsigned NN = ( orderTriang * ( orderTriang + 1)) / 2; // number of Triangulation Points (without merging)

    _ptsTriang = VectorVector( 20 * NN );   // store Triangulation points
    _triangles = VectorVectorU( GetNq() );  // store triangle connection of triangulation points

    SetupPlatonicSolid(); // Get Initial Icosahedron in unit sphere

    for (unsigned id_face = 0; id_face < 20; id_face++){  // loop over every face to triangulate face
        // Triangulate Face -> get triangulation points and triangles (connectors) per face
        Vector a = _vertices[_faces[ id_face ][0]];
        Vector b = _vertices[_faces[ id_face ][1]];
        Vector c = _vertices[_faces[ id_face ][2]];

        unsigned countTP = id_face * NN;
        unsigned countMP = id_face * GetNq() / 20;  // counter for points and connections

        VectorVector p0p1 = Interpolate ( a, b, orderTriang );    // interpolation left side
        VectorVector p0p2 = Interpolate ( a, c, orderTriang );    // interpolation right side

        for (unsigned id = 0; id < orderTriang; id++ )
        {
            VectorVector pts = Interpolate ( p0p1[ id ], p0p2[ id ], id + 1); // interpolate between left and right point in row id

            for (unsigned j = 0 ; j < id + 1 ; j++){
                _ptsTriang[ countTP + j ] = pts[ j ];  // determine Triangulation Points
                // define triangles
                if (id > 0 && j > 0 ){
                    // upwards triangles
                    _triangles[ countMP ] = { countTP + j - 1 - id , countTP + j - 1 , countTP + j };
                    countMP += 1;
                    if (j != id){
                        // downwards triangles
                        _triangles[ countMP ] = { countTP + j - 1 - id , countTP + j - id ,  countTP + j };
                        countMP += 1;
                    }
                }
            }
            countTP += id + 1;
        }
    }
}

void QIcosahedron::SetPointsAndWeights(){

    _points    = VectorVector( GetNq() );   // store Quadrature points (as midpoints of triangles)
    _pointsSphere = VectorVector( GetNq(), Vector( 2 ));
    _weights   = Vector( GetNq() );         // store weights of Quadrature points ( as area of triangles)

    for ( unsigned id = 0; id < GetNq(); id++){ // points are midpoints of subtriangles
        _points[ id ] = 1.0 / 3.0 * ( _ptsTriang[ _triangles[ id ][ 0 ] ] +
                                      _ptsTriang[ _triangles[ id ][ 1 ] ] +
                                      _ptsTriang[ _triangles[ id ][ 2 ] ] );
        _points[ id ] /= norm( _points[ id ] );

        _pointsSphere[ id ][ 0 ] = _points[ id ][ 2 ];  // my = Omega_z
        _pointsSphere[ id ][ 1 ] = atan2( _points[ id ][ 1 ] ,  _points[ id ][ 0 ] );  // phi = acos( Omega_x )

        _weights[ id ] = GetArea( _ptsTriang[ _triangles[ id ][ 0 ] ] ,
                                  _ptsTriang[ _triangles[ id ][ 1 ] ] ,
                                  _ptsTriang[ _triangles[ id ][ 2 ] ] );
    }
}


void QIcosahedron::SetConnectivity() {    // TODO
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

void QIcosahedron::SetNeighbourConnectivity() {

    _neighbours = VectorVectorU(_nq, VectorU( 3u ));
    _neighboursRefined = VectorVectorU(_nqRefined, VectorU( 3u, 0 ));

    // this is not nice programming !!!
    // idea: three nearest points have to be the neighbours

    for (unsigned id = 0; id < _nq; id ++){

        Vector distance1 = Vector( _nq, 0.0 );
        VectorU ids = VectorU( _nq );
        for ( unsigned jd = 0; jd < _nq; jd++ ){
            if ( jd == id ) distance1[ jd ] = 10e5;
            else
                distance1[ jd ] = blaze::l2Norm( _points[id] - _points[jd] );
        }

        std::iota( ids.begin(), ids.end(), 0);
        std::sort( ids.begin(), ids.end(), [&](int i,int j){return distance1[i] < distance1[j];} );
        _neighbours[id] = {ids[0], ids[1], ids[2]};

    }


    /* nott needed yet
    if ( _settings->GetIsLocalRefine() ){
        for (unsigned id = 0; id < _nqRefined; id ++){

            Vector distance1 = Vector( _nqRefined, 0.0 );
            VectorU ids = VectorU( _nqRefined );
            for ( unsigned jd = 0; jd < _nqRefined; jd++ ){
                if ( jd == id ) distance1[ jd ] = 10e5;
                else
                    distance1[ id ] = blaze::l2Norm( _pointsRefined[id] - _pointsRefined[jd] );
            }

            std::iota( ids.begin(), ids.end(), 0);
            std::sort( ids.begin(), ids.end(), [&](int i,int j){return distance1[i] < distance1[j];} );

            _neighboursRefined[id] = {ids[0], ids[1], ids[2]};
        }
    }
    */
}

void QIcosahedron::Refine(){

    _pointsRefined  = VectorVector( GetNqRefined() );
    _weightsRefined = Vector( GetNqRefined() );
    _refineVector   = VectorVectorU( GetNqRefined() );


    for (unsigned id = 0; id < GetNq(); id++){

        Vector a = _ptsTriang[_triangles[ id ][0]];
        Vector b = _ptsTriang[_triangles[ id ][1]];
        Vector c = _ptsTriang[_triangles[ id ][2]];

        if (_settings->GetQuadOrderFine() == _order + 2){
            VectorVectorU triang( 16, VectorU( 3 ) );
            VectorVector p0p1 = Interpolate ( a, b, 5 );    // interpolation left side
            VectorVector p0p2 = Interpolate ( a, c, 5 );    // interpolation right side
            VectorVector p1p2 = Interpolate ( b, c, 5 );    // interpolation down side

            VectorVector ptsTri(15, Vector(3));
            unsigned cT = 0;
            unsigned cc = 0;

            for (unsigned id = 0; id < 5; id++ )
            {
                VectorVector pts = Interpolate ( p0p1[ id ], p0p2[ id ], id + 1); // interpolate between left and right point in row id
                for (unsigned j = 0 ; j < id + 1 ; j++){
                    ptsTri[ cT + j ] = pts[ j ];  // determine Triangulation Points
                    // define triangles
                    if (id > 0 && j > 0 ){
                        triang[ cc ] = { cT + j - 1 - id , cT + j - 1 , cT + j }; // upwards triangles
                        cc += 1;
                        if (j != id){
                            triang[ cc ] = { cT + j - 1 - id , cT + j - id ,  cT + j }; // downwards triangles
                            cc += 1;
                        }
                    }
                }
                cT += id + 1;
            }
            _refineVector[ id ] = Vector( 16u );
            for ( unsigned i = 0; i < 16; i++){
                _pointsRefined[ 16 * id + i ] = 1.0 / 3.0 * ( ptsTri[triang[i][0]] + ptsTri[triang[i][1]] + ptsTri[triang[i][2]] );
                _weightsRefined[ 16 * id + i ] = GetArea ( ptsTri[triang[i][0]], ptsTri[triang[i][1]], ptsTri[triang[i][2]] );
                _refineVector[ id ][ i ] = 16 * id + i;
            }
        }
        else if ( _settings->GetQuadOrderFine() == _order + 1 ){
            VectorVectorU triang( 4, VectorU( 3 ) );
            VectorVector p0p1 = Interpolate ( a, b, 3 );    // interpolation left side
            VectorVector p0p2 = Interpolate ( a, c, 3 );    // interpolation right side
            VectorVector p1p2 = Interpolate ( b, c, 3 );    // interpolation down side

            _pointsRefined[ 4 * id ] = 1.0 / 3.0 * ( p0p1[1] + p0p2[1] + p1p2[1] );         // old quadpoint
            _pointsRefined[ 4 * id + 1 ] = 1.0 / 3.0 * ( p0p1[0] + p0p1[1] + p0p2[1] );     // upper triangle
            _pointsRefined[ 4 * id + 2 ] = 1.0 / 3.0 * ( p0p1[1] + p1p2[0] + p1p2[1] );     // left triangle
            _pointsRefined[ 4 * id + 3 ] = 1.0 / 3.0 * ( p0p2[1] + p1p2[1] + p1p2[2] );     // right triangle

            _weightsRefined[ 4 * id ] = GetArea ( p0p1[1] , p0p2[1] , p1p2[1] );
            _weightsRefined[ 4 * id + 1 ] = GetArea ( p0p1[0] , p0p1[1] , p0p2[1] );
            _weightsRefined[ 4 * id + 2 ] = GetArea ( p0p1[1] , p1p2[0] , p1p2[1] );
            _weightsRefined[ 4 * id + 3 ] = GetArea ( p0p2[1] , p1p2[1] , p1p2[2] );

            _refineVector[ id ] = { 4 * id, 4 * id + 1, 4 * id + 2, 4 * id + 3 };

        }
        else if( _settings->GetQuadOrderFine() == _order){
            _pointsRefined = _points;
            _weightsRefined = _weights;
            _refineVector = VectorVectorU(_nq, VectorU(1u));
            for (unsigned i = 0; i< _nq; i++)
                _refineVector[i] = i;
        }
    }
    for ( unsigned id = 0; id < GetNqRefined(); id++)
        _pointsRefined[ id ] /= norm( _pointsRefined[ id ]);
}

//------------ Helper ---------------//

void QIcosahedron::SetupPlatonicSolid(){
    double r = (1 + std::sqrt(5)) / 2.0;
    _vertices = {   {0.0,  1.0,  r},  { 0.0, 1.0, -r},  {0.0, -1.0,  r},
                    {0.0, -1.0, -r},  { 1.0,  r, 0.0},  {1.0, -r,  0.0},
                    {-1.0,  r, 0.0},  {-1.0, -r, 0.0},  {r,  0.0,  1.0},
                    { r, 0.0, -1.0},  {-r , 0.0 ,1.0},  {-r ,0.0 ,-1.0} };

    for (unsigned id = 0; id < 12; id++){ // normalize vertices to get unit sphere
        _vertices[id] /= norm(_vertices[id]);
    }

    _faces = {  {0, 2, 8},   {0, 2, 10},  {0, 4, 6},  {0, 6, 10},
                {1, 4, 6},   {1, 6, 11},  {1, 4, 9},  {1, 3, 9},
                {1, 3, 11},  {7, 10, 11}, {3, 7, 11}, {3, 5, 7},
                {3, 5, 9},   {5, 8, 9},   {4, 8, 9},  {0, 4, 8},
                {6, 10, 11}, {2, 5, 8},   {2, 5, 7},  {2, 7, 10} };
}

VectorVector QIcosahedron::Interpolate( Vector a, Vector b, unsigned n ){
    if ( n == 1 )
        return VectorVector( 1, a );

    VectorVector out( n, Vector( 3, 0.0 ) );
    if (n == 2) { out[0] = a; out[1] = b; return out; } // exception bc linspace( 2, 0, 1 ) results in [1,1]+

    blaze::DynamicVector<double, true> aT = trans(a);
    double omega = acos (aT * b );
    Vector t = blaze::linspace ( n, 0.0, 1.0 );

    for (unsigned id = 0; id < n; id++)
    {
        out[id] = a * ( sin( (1.0 - t[id]) * omega ) ) / sin ( omega ) + b * ( sin( ( t[id] ) * omega ) ) / sin ( omega );
    }
    return out;
}

double QIcosahedron::GetArea( Vector a, Vector b, Vector c){

    a /= norm(a); b /= norm(b); c /= norm(c);

    blaze::DynamicVector<double, true> aT = trans(a);
    blaze::DynamicVector<double, true> bT = trans(a);
    blaze::DynamicVector<double, true> cT = trans(a);
    // calc Distances
    double A = acos((bT * c));
    double B = acos((aT * c));
    double C = acos((aT * b));
    // calc Angles
    double alpha = acos((cos(A) - cos(B) * cos(C) ) / (sin(B) * sin(C)));
    double beta =  acos((cos(B) - cos(C) * cos(A) ) / (sin(C) * sin(A)));
    double gamma = acos((cos(C) - cos(A) * cos(B) ) / (sin(A) * sin(B)));

    return (alpha + beta + gamma - PI);
}

void QIcosahedron::ResetValues(){
    _nqRefined = _nq;
    _pointsRefined = _points;
    _weightsRefined = _weights;
    for( unsigned id = 0; id < GetNq(); id++)
        _refineVector[id] = {id};
}

bool QIcosahedron::CheckOrder()
{
    return true;
}