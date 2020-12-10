/*!
 * @file sphericalmonomials.cpp
 * @brief Class for efficient computation of a spherical monomial basis
 * @author S. Schotthöfer
 */

#include "toolboxes/sphericalmonomials.h"

SphericalMonomials::SphericalMonomials( unsigned L_degree ) {
    _LMaxDegree = L_degree;
    _spatialDim = 3;    // Spatial dimension 2 is hardcoded

    _YBasis = Vector( GetBasisSize() );
}

Vector SphericalMonomials::ComputeSphericalBasis( double my, double phi ) {
    unsigned idx_vector = 0;
    // go over all degrees of polynomials
    for( unsigned idx_degree = 0; idx_degree <= _LMaxDegree; idx_degree++ ) {
        // elem = Omega_x^a+ Omega_y^b +Omega_z^c : a+b+c = idx_degree
        double omega_X = Omega_x( my, phi );
        double omega_Y = Omega_y( my, phi );
        double omega_Z = Omega_z( my );
        for( unsigned a = 0; a <= idx_degree; a++ ) {
            for( unsigned b = 0; b <= idx_degree - a; b++ ) {
                unsigned c          = idx_degree - a - b;    // c uniquely defined
                _YBasis[idx_vector] = Power( omega_X, a ) * Power( omega_Y, b ) * Power( omega_Z, c );
                idx_vector++;
            }
        }
    }
    return _YBasis;
}

Vector SphericalMonomials::ComputeSphericalBasis( double x, double y, double z ) {
    // transform (x,y,z) into (my,phi)
    double my  = z;
    double phi = 0.0;

    if( y >= 0 )
        phi = acos( x );
    else
        phi = 2 * M_PI - acos( x );

    return ComputeSphericalBasis( my, phi );
}

unsigned SphericalMonomials::GetBasisSize() {
    unsigned basisLen = 0;
    for( unsigned idx_degree = 0; idx_degree <= _LMaxDegree; idx_degree++ ) {
        basisLen += ComputeDimensionSize( idx_degree );
    }
    return basisLen;
}

unsigned SphericalMonomials::ComputeDimensionSize( unsigned degree ) {
    return Factorial( degree + _spatialDim - 1 ) / ( Factorial( degree ) * Factorial( _spatialDim - 1 ) );
}

unsigned SphericalMonomials::Factorial( unsigned n ) { return ( n == 1 || n == 0 ) ? 1 : Factorial( n - 1 ) * n; }

double SphericalMonomials::Omega_x( double my, double phi ) { return sqrt( 1 - my * my ) * sin( phi ); }
double SphericalMonomials::Omega_y( double my, double phi ) { return sqrt( 1 - my * my ) * cos( phi ); }
double SphericalMonomials::Omega_z( double my ) { return my; }

double SphericalMonomials::Power( double basis, unsigned exponent ) {
    if( exponent == 0 ) return 1.0;
    double result = basis;
    for( unsigned i = 1; i < exponent; i++ ) {
        result = result * basis;
    }
    return result;
}
