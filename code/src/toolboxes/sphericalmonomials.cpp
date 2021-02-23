/*!
 * @file sphericalmonomials.cpp
 * @brief Class for efficient computation of a spherical monomial basis
 * @author S. Schotthöfer
 */

#include "toolboxes/sphericalmonomials.h"
#include "toolboxes/errormessages.h"

SphericalMonomials::SphericalMonomials( unsigned L_degree ) {
    _LMaxDegree = L_degree;
    _spatialDim = 3;    // Spatial dimension 2 is hardcoded

    _YBasis = Vector( GetBasisSize() );
}

SphericalMonomials::SphericalMonomials( unsigned L_degree, unsigned short spatialDim ) {
    _LMaxDegree = L_degree;
    _spatialDim = (unsigned)spatialDim;

    _YBasis = Vector( GetBasisSize() );
}

Vector SphericalMonomials::ComputeSphericalBasis( double my, double phi ) {

    switch( _spatialDim ) {
        case 1: return ComputeSphericalBasis1D( my, phi ); break;
        case 2: return ComputeSphericalBasis2D( my, phi ); break;
        default: return ComputeSphericalBasis3D( my, phi ); break;
    }
}

Vector SphericalMonomials::ComputeSphericalBasis1D( double my, double phi ) {
    unsigned idx_vector = 0;
    unsigned a;
    double omega_X;
    // go over all degrees of polynomials
    for( unsigned idx_degree = 0; idx_degree <= _LMaxDegree; idx_degree++ ) {
        // elem = Omega_x^a : a = idx_degree
        omega_X             = Omega_x( my, phi );
        a                   = idx_degree;    // a uniquely defined
        _YBasis[idx_vector] = Power( omega_X, a );
        idx_vector++;
    }
    return _YBasis;
}

Vector SphericalMonomials::ComputeSphericalBasis2D( double my, double phi ) {
    unsigned idx_vector = 0;
    double omega_X, omega_Y;
    unsigned a, b;

    // go over all degrees of polynomials
    for( unsigned idx_degree = 0; idx_degree <= _LMaxDegree; idx_degree++ ) {
        // elem = Omega_x^a+ Omega_y^b  : a+b = idx_degree
        omega_X = Omega_x( my, phi );
        omega_Y = Omega_y( my, phi );
        for( a = 0; a <= idx_degree; a++ ) {
            b                   = idx_degree - a;    // b uniquely defined
            _YBasis[idx_vector] = Power( omega_X, a ) * Power( omega_Y, b );
            idx_vector++;
        }
    }
    return _YBasis;
}

Vector SphericalMonomials::ComputeSphericalBasis3D( double my, double phi ) {
    unsigned idx_vector = 0;
    unsigned a, b, c;

    double omega_X = Omega_x( my, phi );
    double omega_Y = Omega_y( my, phi );
    double omega_Z = Omega_z( my );
    // go over all degrees of polynomials
    for( unsigned idx_degree = 0; idx_degree <= _LMaxDegree; idx_degree++ ) {
        // elem = Omega_x^a+ Omega_y^b +Omega_z^c : a+b+c = idx_degree
        for( a = 0; a <= idx_degree; a++ ) {
            for( b = 0; b <= idx_degree - a; b++ ) {
                c                   = idx_degree - a - b;    // c uniquely defined
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
        basisLen += GetCurrDegreeSize( idx_degree );
    }
    return basisLen;
}

unsigned SphericalMonomials::GetCurrDegreeSize( unsigned currDegreeL ) {
    return Factorial( currDegreeL + _spatialDim - 1 ) / ( Factorial( currDegreeL ) * Factorial( _spatialDim - 1 ) );
}

unsigned SphericalMonomials::GetGlobalIndexBasis( int l_degree, int k_order ) {
    if( l_degree < 0 ) ErrorMessages::Error( "Negative polynomial degrees not supported.", CURRENT_FUNCTION );
    if( k_order < 0 ) ErrorMessages::Error( "Order k of spherical monomial basis must not be negative.", CURRENT_FUNCTION );
    if( k_order >= (int)GetCurrDegreeSize( l_degree ) )
        ErrorMessages::Error( "Order k of spherical monomial basis out of bounds.", CURRENT_FUNCTION );

    unsigned basisLen = 0;
    for( unsigned idx_degree = 0; idx_degree < (unsigned)l_degree; idx_degree++ ) {
        basisLen += GetCurrDegreeSize( idx_degree );
    }
    return basisLen + (unsigned)k_order;
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
