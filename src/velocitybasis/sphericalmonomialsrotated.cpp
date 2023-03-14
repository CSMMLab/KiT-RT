/*!
 * @file sphericalmonomials.cpp
 * @brief Class for efficient computation of a spherical monomial basis with rotated velocity reference frame
 * @author S. Schotthöfer
 */

#include "velocitybasis/sphericalmonomialsrotated.hpp"
#include "toolboxes/errormessages.hpp"

SphericalMonomialsRotated::SphericalMonomialsRotated( unsigned L_degree ) {
    _LMaxDegree = L_degree;

    _spatialDim = 3;    // Spatial dimension 2 is hardcoded

    _YBasis = Vector( GetBasisSize() );
}

SphericalMonomialsRotated::SphericalMonomialsRotated( unsigned L_degree, unsigned short spatialDim ) {
    _LMaxDegree = L_degree;
    _spatialDim = (unsigned)spatialDim;

    if( _spatialDim != 2 ) ErrorMessages::Error( "Spatial dimension other than 2 not supported for rotated basis." , CURRENT_FUNCTION );
    if( _LMaxDegree > 2 ) ErrorMessages::Error( "Basis degree higher than 2 not supported for rotated basis." , CURRENT_FUNCTION );

    _YBasis = Vector( GetBasisSize() );
}

Vector SphericalMonomialsRotated::ComputeSphericalBasis( double my, double phi, double r ) {
    // default for (maximal) radius r = 1
    switch( _spatialDim ) {
        case 1: return ComputeSphericalBasis1D( my, r ); break;
        case 2: return ComputeSphericalBasis2D( my, phi, r ); break;
        default: return ComputeSphericalBasis3D( my, phi, r ); break;
    }
}

Vector SphericalMonomialsRotated::ComputeSphericalBasis1D( double my, double r ) {
    // default for (maximal) radius r = 1
    unsigned idx_vector = 0;
    unsigned a;
    double omega_Z;
    // go over all degrees of polynomials
    for( unsigned idx_degree = 0; idx_degree <= _LMaxDegree; idx_degree++ ) {
        // elem = Omega_x^a : a = idx_degree
        omega_Z             = Omega_z( my, r );
        a                   = idx_degree;    // a uniquely defined
        _YBasis[idx_vector] = Power( omega_Z, a );
        idx_vector++;
    }
    return _YBasis;
}

Vector SphericalMonomialsRotated::ComputeSphericalBasis2D( double my, double phi, double r ) {
    // default for (maximal) radius r = 1
    unsigned idx_vector = 0;
    double omegaX, omegaY;
    unsigned a, b;

    omegaX = Omega_x( my, phi, r );
    omegaY = Omega_y( my, phi, r );

    Vector omegaXY{ omegaX, omegaY };
    Matrix rotationMatrix = CreateRotator( omegaXY );
    Matrix rotationMatrixTrans  = blaze::trans( rotationMatrix );

    omegaXY = RotateM1( omegaXY, rotationMatrix );    // Rotate velocity frame to x-axis

    omegaX = omegaXY[0];
    omegaY = omegaXY[1];

    // go over all degrees of polynomials
    for( unsigned idx_degree = 0; idx_degree <= _LMaxDegree; idx_degree++ ) {
        // elem = Omega_x^a+ Omega_y^b  : a+b = idx_degree

        for( a = 0; a <= idx_degree; a++ ) {
            b                   = idx_degree - a;    // b uniquely defined
            _YBasis[idx_vector] = Power( omegaX, a ) * Power( omegaY, b );
            idx_vector++;
        }
    }

    Vector m1{ _YBasis[1], _YBasis[2] };
    Matrix m2{ { _YBasis[3], 0.5 * _YBasis[4] }, { 0.5 * _YBasis[4], _YBasis[5] } };

    // Rotate basis back
    m1 = RotateM1( m1, rotationMatrixTrans );
    m2 = RotateM2( m2, rotationMatrixTrans, rotationMatrix );

    _YBasis[1] = m1[0];
    _YBasis[2] = m1[1];
    _YBasis[3] = m2( 0, 0 );
    _YBasis[4] = 2 * m2( 1, 0 );
    _YBasis[5] = m2( 1, 1 );

    //for( unsigned idx_degree = 0; idx_degree <= _LMaxDegree; idx_degree++ ) {
    //    // elem = Omega_x^a+ Omega_y^b  : a+b = idx_degree
    //    omegaX = Omega_x( my, phi, r );
    //    omegaY = Omega_y( my, phi, r );
    //    for( a = 0; a <= idx_degree; a++ ) {
    //        b                   = idx_degree - a;    // b uniquely defined
    //        _YBasis[idx_vector] = Power( omegaX, a ) * Power( omegaY, b );
    //        idx_vector++;
    //    }
    //}

    return _YBasis;
}

Vector SphericalMonomialsRotated::ComputeSphericalBasis3D( double my, double phi, double r ) {
    // default for radius r = 1
    unsigned idx_vector = 0;
    unsigned a, b, c;

    double omega_X = Omega_x( my, phi, r );
    double omega_Y = Omega_y( my, phi, r );
    double omega_Z = Omega_z( my, r );
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

Vector SphericalMonomialsRotated::ComputeSphericalBasisKarthesian( double x, double y, double z ) {
    // transform (x,y,z) into (my,phi)
    double my  = z;                                // my = z
    double phi = atan2( y, x );                    // phi in [-pi,pi]
    double r   = sqrt( x * x + y * y + z * z );    // radius r

    // adapt intervall s.t. phi in [0,2pi]
    if( phi < 0 ) {
        phi = 2 * M_PI + phi;
    }

    return ComputeSphericalBasis( my, phi, r );
}

unsigned SphericalMonomialsRotated::GetBasisSize() {
    unsigned basisLen = 0;
    for( unsigned idx_degree = 0; idx_degree <= _LMaxDegree; idx_degree++ ) {
        basisLen += GetCurrDegreeSize( idx_degree );
    }
    return basisLen;
}

unsigned SphericalMonomialsRotated::GetCurrDegreeSize( unsigned currDegreeL ) {
    return Factorial( currDegreeL + _spatialDim - 1 ) / ( Factorial( currDegreeL ) * Factorial( _spatialDim - 1 ) );
}

unsigned SphericalMonomialsRotated::GetGlobalIndexBasis( int l_degree, int k_order ) {
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

unsigned SphericalMonomialsRotated::Factorial( unsigned n ) { return ( n == 1 || n == 0 ) ? 1 : Factorial( n - 1 ) * n; }

double SphericalMonomialsRotated::Omega_x( double my, double phi, double r ) { return r * sqrt( 1 - my * my ) * cos( phi ); }
double SphericalMonomialsRotated::Omega_y( double my, double phi, double r ) { return r * sqrt( 1 - my * my ) * sin( phi ); }
double SphericalMonomialsRotated::Omega_z( double my, double r ) { return r * my; }

double SphericalMonomialsRotated::Power( double basis, unsigned exponent ) {
    if( exponent == 0 ) return 1.0;               // basis^0
    double result = basis;                        // basis^1
    for( unsigned i = 1; i < exponent; i++ ) {    // exp> 1
        result = result * basis;
    }
    return result;
}

Vector SphericalMonomialsRotated::RotateM1( Vector& vec, Matrix& R ) { return R * vec; }

Matrix SphericalMonomialsRotated::RotateM2( Matrix& vec, Matrix& R, Matrix& Rt ) { return R * vec * Rt; }

Matrix SphericalMonomialsRotated::CreateRotator( const Vector& uFirstMoment ) {
    double a = uFirstMoment[0];
    double b = uFirstMoment[1];
    double c, s, r;

    r = norm( uFirstMoment );    // sqrt( a * a + b * b );
    c = a / r;
    s = -b / r;

    return Matrix{ { c, -s }, { s, c } };    // Rotation Matrix
}
