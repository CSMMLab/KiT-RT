//
// Created by chinsp on 27/10/21.
//

#include "../../include/quadratures/qicosahedrontriang.h"
#include "common/config.h"
#include "toolboxes/errormessages.h"

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
Uses triangulation vertex points as Quadrature Points
Uses sum of areas of corresponding subtriangles as Quadrature Weights
Refinement not defined
*/


QIcosahedronII::QIcosahedronII( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    CheckOrder();
    SetNq();
    SetupTriangulation();
    SetPointsAndWeights();
    SetConnectivity();
}

QIcosahedronII::QIcosahedronII( unsigned order ) : QuadratureBase( order ) {
    SetName();
    CheckOrder();
    SetNq();
    SetupTriangulation();
    SetPointsAndWeights();
    SetConnectivity();
}

void QIcosahedronII::SetNq() {
    unsigned orderTriang = pow( 2, _order-1) + 1;
    _nq =  20 * ( orderTriang * ( orderTriang + 1)) / 2; // without merging
    _nq -= 4 * 12; // every vertex point is used in 5 faces
    _nq -= 20 * ( orderTriang - 2 ) * 3 / 2 ; // every edge point is used in two faces
}


void QIcosahedronII::SetupTriangulation() {

    unsigned orderTriang = pow( 2, _order-1) + 1;
    unsigned NN = ( orderTriang * ( orderTriang + 1)) / 2; // number of Triangulation Points (without merging)

    _ptsTriang = VectorVector( 20 * NN );   // store Triangulation points
    _triangles = VectorVectorU( 20 * pow(4, _order -1), Vector( 3u ) );  // store triangle connection of triangulation points

    SetupPlatonicSolid( ); // Get Initial Icosahedron in unit sphere


    for (unsigned id_face = 0; id_face < _faces.size(); id_face++){  // loop over every face to triangulate face
        // Triangulate Face -> get triangulation points and triangles (connectors) per face
        Vector a = _vertices[_faces[ id_face ][0]];
        Vector b = _vertices[_faces[ id_face ][1]];
        Vector c = _vertices[_faces[ id_face ][2]];

        unsigned countTP = id_face * NN;
        unsigned countMP = id_face * pow(4, _order - 1);  // counter for points and connections

        VectorVector p0p1 = Interpolate ( a, b, orderTriang );    // interpolation left side
        VectorVector p0p2 = Interpolate ( a, c, orderTriang );    // interpolation right side

        for (unsigned id = 0; id < orderTriang; id++ )
        {
            VectorVector pts = Interpolate ( p0p1[ id ], p0p2[ id ], id + 1 ); // interpolate between left and right point in row id

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

void QIcosahedronII::SetPointsAndWeights(){

    _points    = VectorVector( GetNq() );   // store Quadrature points (as merged triangulation points)
    _pointsSphere = VectorVector( GetNq() , Vector( 2u ));
    _weights   = Vector( GetNq() );         // store weights of Quadrature points ( as area of connected triangles)
    VectorVectorU matches = VectorVectorU( _ptsTriang.size() * _ptsTriang.size(), Vector (2u));
    double wsum = 0;


    // merge duplicate points
    double ind = 0;
    for ( unsigned id = 0; id < _ptsTriang.size(); id ++){
        for ( unsigned jd = 0; jd < id ; jd ++){
            // jd == id not reached :)
            if ( norm( _ptsTriang[id] - _ptsTriang[jd] ) < 1e-10 ){
                matches[ ind ] = { jd, id } ; // smaller value is in front
                ind += 1;
            }
        }
    }

    matches.resize( ind );
    for ( unsigned id = matches.size() - 1 ; id < matches.size(); id--){ // !!! (unsigned) -1 is \infty
        _ptsTriang[ matches[id][1] ] = {nan(""), nan(""), nan("")} ;
        // adapt triangles
        for ( unsigned jd = 0; jd < _triangles.size(); jd++ ){
            for ( unsigned kd = 0; kd < _triangles[jd].size(); kd++ ){
                if ( _triangles[ jd ][ kd ] == matches[ id ][ 1 ]){
                    _triangles[ jd ][ kd ] = matches[ id ][ 0 ] ;
                }
            }
        }
    }

    // Get unique quadPoints and calculate weights
    ind = 0;
    for ( unsigned id = 0; id < _ptsTriang.size(); id++){

        if ( blaze::isnan( _ptsTriang[ id ][0] ) ) { continue; }
        else {
            _points[ ind ] = _ptsTriang[ id ];

            _pointsSphere[ ind ][ 0 ] = _points[ ind ][ 2 ];  // my = Omega_z
            _pointsSphere[ ind ][ 1 ] = atan2( _points[ ind ][ 1 ] ,  _points[ ind ][ 0 ] );  // phi = acos( Omega_x )

            for ( unsigned jd = 0; jd < _triangles.size(); jd++ ){
                for ( unsigned kd = 0; kd < _triangles[0].size(); kd++ ){
                    if ( _triangles[ jd ][ kd ] == id ){
                        // Calculate area of the triangle containing point id
                        double area = GetArea( _ptsTriang[ _triangles[ jd ][ 0 ] ], _ptsTriang[ _triangles[ jd ][ 1 ] ], _ptsTriang[ _triangles[ jd ][ 2 ] ] );
                        // add up to weights
                        _weights[ ind ] += area;
                        break;
                    }
                }
            }
            // wsum += _weights[ ind ]; // later used to normalize
            ind += 1;
        }
    }
    // Normalize Weights to get sum( w(x) ) = 1.
    for ( unsigned id = 0; id < _weights.size(); id ++ )
        _weights[ id ] /= 3.0; // !!!!!!!!!!!

}


void QIcosahedronII::SetConnectivity() {    // TODO
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

//------------ Helper ---------------//

void QIcosahedronII::SetupPlatonicSolid(){
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

VectorVector QIcosahedronII::Interpolate( Vector a, Vector b, unsigned n ){
    if ( n == 1 )
        return VectorVector( 1, a );

    VectorVector out( n, Vector( 3, 0.0 ) );
    if (n == 2) { out[0] = a; out[1] = b; return out; } // exception bc linspace( 2, 0, 1 ) results in [1,1]

    blaze::DynamicVector<double, true> aT = trans(a);
    double omega = acos (aT * b );
    Vector t = blaze::linspace ( n, 0.0, 1.0 );

    for (unsigned id = 0; id < n; id++)
    {
        out[id] = a * ( sin( (1.0 - t[id]) * omega ) ) / sin ( omega ) + b * ( sin( ( t[id] ) * omega ) ) / sin ( omega );
    }
    return out;
}

double QIcosahedronII::GetArea( Vector a, Vector b, Vector c ){

    // weights are surface areas of spherical triangles
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



bool QIcosahedronII::CheckOrder()
{
    return true;
}