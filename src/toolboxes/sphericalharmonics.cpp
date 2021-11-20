#include "toolboxes/sphericalharmonics.hpp"
#include "toolboxes/errormessages.hpp"
#include <math.h>

SphericalHarmonics::SphericalHarmonics( unsigned L_degree ) {
    _LMaxDegree = L_degree;

    unsigned associatedLegendrePolySize = GlobalIdxAssLegendreP( _LMaxDegree, _LMaxDegree ) + 1;

    _aParam       = std::vector<double>( associatedLegendrePolySize, 0.0 );
    _bParam       = std::vector<double>( associatedLegendrePolySize, 0.0 );
    _assLegendreP = std::vector<double>( associatedLegendrePolySize, 0.0 );

    ComputeCoefficients();

    _YBasis = Vector( GetBasisSize(), 0.0 );
}

unsigned SphericalHarmonics::GetBasisSize() { return GetGlobalIndexBasis( _LMaxDegree, _LMaxDegree ) + 1; /* +1, since globalIdx computes indices */ }

unsigned SphericalHarmonics::GetCurrDegreeSize( unsigned currDegreeL ) { return 2 * currDegreeL + 1; }

unsigned SphericalHarmonics::GetGlobalIndexBasis( int l_degree, int k_order ) {
    if( l_degree < 0 ) ErrorMessages::Error( "Negative polynomial degrees not supported.", CURRENT_FUNCTION );
    if( k_order < -l_degree || k_order > l_degree )
        ErrorMessages::Error( "Order k of spherical harmonics basis must be in bounds -l<= k <= l.", CURRENT_FUNCTION );

    return (unsigned)( k_order + l_degree /*number of previous indices in current level*/ +
                       l_degree * l_degree ) /* number of previous indices untill level l-1 */;
}

Vector SphericalHarmonics::ComputeSphericalBasis( double my, double phi, double r ) {
    // radius r is not used and always assumed to be 1

    ComputeAssLegendrePoly( my );
    ComputeYBasis( phi );
    return _YBasis;
}

Vector SphericalHarmonics::ComputeSphericalBasisKarthesian( double x, double y, double z ) {

    // transform (x,y,z) into (my,phi)
    double my  = z;
    double phi = 0.0;

    if( y >= 0 )
        phi = acos( x );
    else
        phi = 2 * M_PI - acos( x );

    ComputeAssLegendrePoly( my );
    ComputeYBasis( phi );
    return _YBasis;
}

std::vector<double> SphericalHarmonics::GetAssLegendrePoly( const double my ) {
    ComputeAssLegendrePoly( my );
    return _assLegendreP;
}

void SphericalHarmonics::ComputeCoefficients() {
    // m in paper is here denoted by k
    double ls   = 0.0;    // l^2
    double lm1s = 0.0;    // (l-1)^2
    double ks   = 0.0;    // k^2
    for( unsigned l_idx = 2; l_idx <= _LMaxDegree; l_idx++ ) {
        ls   = l_idx * l_idx;
        lm1s = ( l_idx - 1 ) * ( l_idx - 1 );
        for( unsigned k_idx = 0; k_idx < l_idx - 1; k_idx++ ) {
            ks = k_idx * k_idx;

            _aParam[GlobalIdxAssLegendreP( l_idx, k_idx )] = std::sqrt( ( 4 * ls - 1. ) / ( ls - ks ) );
            _bParam[GlobalIdxAssLegendreP( l_idx, k_idx )] = -std::sqrt( ( lm1s - ks ) / ( 4 * lm1s - 1. ) );
        }
    }
}

void SphericalHarmonics::ComputeAssLegendrePoly( const double my ) {
    const double sintheta = std::sqrt( 1. - my * my );
    double temp           = std::sqrt( .5 / M_PI );

    _assLegendreP[GlobalIdxAssLegendreP( 0, 0 )] = temp;

    if( _LMaxDegree > 0 ) {
        const double SQRT3                           = std::sqrt( 3 );       // 1.7320508075688772935
        _assLegendreP[GlobalIdxAssLegendreP( 1, 0 )] = my * SQRT3 * temp;    // 1.224744871391589
        const double SQRT3DIV2                       = -std::sqrt( 3. / 2. );

        temp                                         = SQRT3DIV2 * sintheta * temp;
        _assLegendreP[GlobalIdxAssLegendreP( 1, 1 )] = temp;

        for( unsigned l_idx = 2; l_idx <= _LMaxDegree; l_idx++ ) {
            for( unsigned k_idx = 0; k_idx < l_idx - 1; k_idx++ ) {
                _assLegendreP[GlobalIdxAssLegendreP( l_idx, k_idx )] =
                    _aParam[GlobalIdxAssLegendreP( l_idx, k_idx )] *
                    ( my * _assLegendreP[GlobalIdxAssLegendreP( l_idx - 1, k_idx )] +
                      _bParam[GlobalIdxAssLegendreP( l_idx, k_idx )] * _assLegendreP[GlobalIdxAssLegendreP( l_idx - 2, k_idx )] );
            }
            _assLegendreP[GlobalIdxAssLegendreP( l_idx, l_idx - 1 )] = my * std::sqrt( 2 * ( l_idx - 1 ) + 3 ) * temp;

            temp = -std::sqrt( 1.0 + 0.5 / l_idx ) * sintheta * temp;

            _assLegendreP[GlobalIdxAssLegendreP( l_idx, l_idx )] = temp;
        }
    }
}

void SphericalHarmonics::ComputeYBasis( const double phi ) {
    for( unsigned l_idx = 0; l_idx <= _LMaxDegree; l_idx++ ) {
        _YBasis[GetGlobalIndexBasis( l_idx, 0 )] = _assLegendreP[GlobalIdxAssLegendreP( l_idx, 0 )] * 0.5 * M_SQRT2;    // M_SQRT2 = sqrt(2)
    }

    // helper constants
    double c1 = 1.0;
    double c2 = std::cos( phi );
    double s1 = 0.0;
    double s2 = -std::sin( phi );
    double tc = 2.0 * c2;

    double s = 0.0;
    double c = 0.0;

    for( unsigned k_idx = 1; k_idx <= _LMaxDegree; k_idx++ ) {
        s  = tc * s1 - s2;    // addition theorem
        c  = tc * c1 - c2;    // addition theorem
        s2 = s1;
        s1 = s;
        c2 = c1;
        c1 = c;
        for( unsigned l_idx = k_idx; l_idx <= _LMaxDegree; l_idx++ ) {
            _YBasis[GetGlobalIndexBasis( l_idx, -k_idx )] = _assLegendreP[GlobalIdxAssLegendreP( l_idx, k_idx )] * s;
            _YBasis[GetGlobalIndexBasis( l_idx, k_idx )]  = _assLegendreP[GlobalIdxAssLegendreP( l_idx, k_idx )] * c;
        }
    }
}
