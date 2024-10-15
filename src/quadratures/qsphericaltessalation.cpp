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

std::vector<std::array<std::array<double, 3>, 3>> QSphericalTessalation::generate_tessellation( const std::array<std::array<double, 3>, 3>& triangle,
                                                                                                int order ) {
    std::vector<std::array<std::array<double, 3>, 3>> triangles;

    // Vertices of the triangle
    std::array<double, 3> A = triangle[0];
    std::array<double, 3> B = triangle[1];
    std::array<double, 3> C = triangle[2];

    // Generate grid points along edges and inside the triangle
    std::vector<std::array<double, 3>> grid_points;
    for( int i = 0; i <= order; ++i ) {
        for( int j = 0; j <= order - i; ++j ) {
            double alpha                = static_cast<double>( i ) / order;
            double beta                 = static_cast<double>( j ) / order;
            double gamma                = 1.0 - alpha - beta;
            std::array<double, 3> point = {
                alpha * A[0] + beta * B[0] + gamma * C[0], alpha * A[1] + beta * B[1] + gamma * C[1], alpha * A[2] + beta * B[2] + gamma * C[2] };
            grid_points.push_back( point );
        }
    }

    // Create triangles from the grid points
    auto idx = [&](int i, int j) -> int {
        if (i < 0 || j < 0 || i > order || j > order - i) {
            throw std::out_of_range("Index out of range in idx calculation");
        }
        return i * (order + 1) - (i * (i - 1)) / 2 + j;
    };
    for( int i = 0; i < order; ++i ) {
        for( int j = 0; j < order - i; ++j ) {
            std::array<double, 3> p1 = grid_points[idx( i, j )];
            std::array<double, 3> p2 = grid_points[idx( i + 1, j )];
            std::array<double, 3> p3 = grid_points[idx( i, j + 1 )];
            triangles.push_back( { p1, p2, p3 } );

            if( j < order - i - 1 ) {
                std::array<double, 3> p4 = grid_points[idx( i + 1, j + 1 )];
                triangles.push_back( { p2, p4, p3 } );
            }
        }
    }
    return triangles;
}

void QSphericalTessalation::SetPointsAndWeights() {
    // Define basal triangle vertices
    std::array<double, 3> A = { 1.0, 0.0, 0.0 };
    std::array<double, 3> B = { 0.0, 1.0, 0.0 };
    std::array<double, 3> C = { 0.0, 0.0, 1.0 };

    // Generate the tessellation for the given order
    std::vector<std::array<std::array<double, 3>, 3>> final_triangles = generate_tessellation( { A, B, C }, _order );

    // Compute centroids of triangles and map them to the unit sphere
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

void QSphericalTessalation::SetNq() { _nq = 4 * pow( GetOrder(), 2 ); }

void QSphericalTessalation::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

bool QSphericalTessalation::CheckOrder() { return true; }

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