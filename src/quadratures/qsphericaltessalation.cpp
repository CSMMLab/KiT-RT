#include "quadratures/qsphericaltessalation.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

QSphericalTessalation::QSphericalTessalation( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    CheckOrder();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
    _supportedDimensions = { 2 };    // Extension to 3d is straightforward
}

void QSphericalTessalation::SetPointsAndWeights() {
    // Define basal triangle vertices
    std::array<double, 3> A                             = { 1.0, 0.0, 0.0 };
    std::array<double, 3> B                             = { 0.0, 1.0, 0.0 };
    std::array<double, 3> C                             = { 0.0, 0.0, 1.0 };
    std::array<std::array<double, 3>, 3> basal_triangle = { A, B, C };

    // Perform recursive subdivision
    std::vector<std::array<std::array<double, 3>, 3>> final_triangles = recursive_subdivide( { basal_triangle }, _order );

    // Compute centroids and map to unit sphere
    std::vector<std::array<double, 3>> centroids;
    std::vector<std::array<double, 3>> mapped_points;
    std::vector<double> weights;

    for( const auto& tri : final_triangles ) {
        auto centroid        = compute_centroid( tri );
        auto mapped_centroid = map_to_unit_sphere( centroid );
        centroids.push_back( centroid );
        mapped_points.push_back( mapped_centroid );

        // Map triangle vertices to unit sphere
        std::array<double, 3> a = map_to_unit_sphere( tri[0] );
        std::array<double, 3> b = map_to_unit_sphere( tri[1] );
        std::array<double, 3> c = map_to_unit_sphere( tri[2] );

        // Calculate area using spherical excess
        double area = spherical_triangle_area( a, b, c );
        weights.push_back( area );
    }

    // Vectors to hold the full set of points and weights for the upper hemisphere
    std::vector<std::array<double, 3>> full_points;
    std::vector<double> full_weights;

    // Perform reflection and permutation
    reflect_and_permute( mapped_points, weights, full_points, full_weights );

    _nq = full_points.size();
    _pointsKarth.resize( _nq );
    _pointsSphere.resize( _nq );
    _weights.resize( _nq );
    for( size_t i = 0; i < _nq; ++i ) {
        _pointsKarth[i] = { full_points[i][0], full_points[i][1], full_points[i][2] };
        _weights[i]     = full_weights[i];
    }
    double w = 0.0;
    // Transform _points to _pointsSphere ==>transform (x,y,z) into (my,phi)
    for( unsigned idx = 0; idx < _nq; idx++ ) {
        _pointsSphere[idx].resize( 3 );                                                 // (my,phi)
        _pointsSphere[idx][0] = _pointsKarth[idx][2];                                   // my = z
        _pointsSphere[idx][1] = atan2( _pointsKarth[idx][1], _pointsKarth[idx][0] );    // phi in [-pi,pi]
        _pointsSphere[idx][2] = 1.0;                                                    // radius r

        // adapt intervall s.t. phi in [0,2pi]
        if( _pointsSphere[idx][1] < 0 ) {
            _pointsSphere[idx][1] = 2 * M_PI + _pointsSphere[idx][1];
        }

        // std::cout << _pointsKarth[idx][0] << " " << _pointsKarth[idx][1] << " " << _pointsKarth[idx][2] << std::endl;
        // std::cout << _pointsSphere[idx][0] << " " << _pointsSphere[idx][1] << " " << _pointsSphere[idx][2] << std::endl;
        // std::cout << _weights[idx] << std::endl;
        w += _weights[idx];
    }
    // std::cout << w << std::endl;
    // exit( 1 );
}

void QSphericalTessalation::SetNq() { _nq = pow( GetOrder(), 4 ); }

void QSphericalTessalation::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

bool QSphericalTessalation::CheckOrder() { return true; }

// Helper function to compute midpoint of two vectors
std::array<double, 3> QSphericalTessalation::midpoint( const std::array<double, 3>& p1, const std::array<double, 3>& p2 ) {
    return { ( p1[0] + p2[0] ) / 2, ( p1[1] + p2[1] ) / 2, ( p1[2] + p2[2] ) / 2 };
}

// Subdivide a triangle into four smaller triangles
std::vector<std::array<std::array<double, 3>, 3>> QSphericalTessalation::subdivide_triangle( const std::array<std::array<double, 3>, 3>& triangle ) {
    std::array<double, 3> A = triangle[0];
    std::array<double, 3> B = triangle[1];
    std::array<double, 3> C = triangle[2];

    // Compute midpoints of edges
    std::array<double, 3> AB_mid = midpoint( A, B );
    std::array<double, 3> BC_mid = midpoint( B, C );
    std::array<double, 3> CA_mid = midpoint( C, A );

    // Define four new triangles
    std::vector<std::array<std::array<double, 3>, 3>> new_triangles = {
        { A, AB_mid, CA_mid }, { AB_mid, B, BC_mid }, { CA_mid, BC_mid, C }, { AB_mid, BC_mid, CA_mid } };

    return new_triangles;
}

// Recursively subdivide triangles to the desired order
std::vector<std::array<std::array<double, 3>, 3>>
QSphericalTessalation::recursive_subdivide( const std::vector<std::array<std::array<double, 3>, 3>>& triangles, int order ) {
    if( order == 1 ) {
        return triangles;
    }
    else {
        std::vector<std::array<std::array<double, 3>, 3>> new_triangles;
        for( const auto& triangle : triangles ) {
            std::vector<std::array<std::array<double, 3>, 3>> subdivided = subdivide_triangle( triangle );
            new_triangles.insert( new_triangles.end(), subdivided.begin(), subdivided.end() );
        }
        return recursive_subdivide( new_triangles, order - 1 );
    }
}

// Compute the centroid of a triangle
std::array<double, 3> QSphericalTessalation::compute_centroid( const std::array<std::array<double, 3>, 3>& triangle ) {
    std::array<double, 3> A = triangle[0];
    std::array<double, 3> B = triangle[1];
    std::array<double, 3> C = triangle[2];

    return { ( A[0] + B[0] + C[0] ) / 3, ( A[1] + B[1] + C[1] ) / 3, ( A[2] + B[2] + C[2] ) / 3 };
}

// Normalize a vector to map it to the unit sphere
std::array<double, 3> QSphericalTessalation::map_to_unit_sphere( const std::array<double, 3>& point ) {
    double norm = std::sqrt( point[0] * point[0] + point[1] * point[1] + point[2] * point[2] );
    return { point[0] / norm, point[1] / norm, point[2] / norm };
}

// Helper function to calculate dot product between two vectors
double QSphericalTessalation::dot_product( const std::array<double, 3>& v1, const std::array<double, 3>& v2 ) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];    // Correct dot product.
}

// Calculate angle at vertex A of spherical triangle ABC
double
QSphericalTessalation::angle_between_vectors( const std::array<double, 3>& a, const std::array<double, 3>& b, const std::array<double, 3>& c ) {
    double ab = dot_product( b, a );
    double ac = dot_product( c, a );
    double bc = dot_product( b, c );
    return std::acos( ( ab - ac * bc ) / ( std::sqrt( 1 - ac * ac ) * std::sqrt( 1 - bc * bc ) ) );
}

// Calculate area of spherical triangle using spherical excess formula
double
QSphericalTessalation::spherical_triangle_area( const std::array<double, 3>& a, const std::array<double, 3>& b, const std::array<double, 3>& c ) {
    double angle_A = angle_between_vectors( b, a, c );
    double angle_B = angle_between_vectors( c, b, a );
    double angle_C = angle_between_vectors( a, c, b );

    // Spherical excess
    return ( angle_A + angle_B + angle_C ) - M_PI;
}

// Function to reflect and permute points from one octant to generate points and weights for the upper half-sphere
void QSphericalTessalation::reflect_and_permute( const std::vector<std::array<double, 3>>& points,
                                                 const std::vector<double>& weights,
                                                 std::vector<std::array<double, 3>>& full_points,
                                                 std::vector<double>& full_weights ) {
    // Loop through the points in the first octant
    for( size_t i = 0; i < points.size(); ++i ) {
        const auto& point = points[i];
        double weight     = weights[i];

        double x = point[0];
        double y = point[1];
        double z = point[2];

        // Generate the four permutations for the upper hemisphere
        std::array<std::array<double, 3>, 4> permutations = { std::array<double, 3>{ x, y, z },
                                                              std::array<double, 3>{ -x, y, z },
                                                              std::array<double, 3>{ x, -y, z },
                                                              std::array<double, 3>{ -x, -y, z } };

        // Append each reflected point and its corresponding weight
        for( const auto& perm : permutations ) {
            full_points.push_back( perm );
            full_weights.push_back( weight );
        }
    }

    // Optional: Project the points onto the z=0 plane (for further use, if necessary)
    // for( auto& p : full_points ) {
    //    p[2] = 0.0;    // Set z to 0 for all points
    //}
}