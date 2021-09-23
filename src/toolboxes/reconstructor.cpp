#include "toolboxes/reconstructor.h"
#include "common/config.h"

Reconstructor::Reconstructor( Config* settings ) { _reconsOrder = settings->GetReconsOrder(); }

double FortSign( double a, double b ) {
    if( b > 0.0 ) return std::fabs( a );
    if( b < 0.0 ) return -std::fabs( a );
    return 0.0;
}

template <typename T> double sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double LMinMod( double sL, double sR ) { return 0.5 * ( FortSign( 1.0, sL ) + FortSign( 1., sR ) ) * fmin( std::fabs( sL ), std::fabs( sR ) ); }

double LMinMod3( double s1, double s2, double s3 ) {
    auto sn1 = sgn(s1);
    auto sn2 = sgn(s2);
    auto sn3 = sgn(s3);

    double s = 0.0;
    if(sn1 == sn2 && sn2 == sn3){
        s = sn1 * fmin(fmin( std::fabs( s1 ), std::fabs( s2 ) ), std::fabs( s3 ) );
    }

    return s;
}

double LVanLeer( double sL, double sR ) {

    return ( FortSign( 1.0, sL ) + FortSign( 1.0, sR ) ) * std::fabs( sL ) * std::fabs( sR ) / ( std::fabs( sL ) + std::fabs( sR ) + 0.0000001 );
}

double LSuperBee( double sL, double sR ) {

    if( sR >= 0.5 * sL && sR <= 2.0 * sL ) return 0.5 * ( FortSign( 1.0, sL ) + FortSign( 1., sR ) ) * fmax( std::fabs( sL ), std::fabs( sR ) );
    if( sR < 0.5 * sL && sR > 2.0 * sL ) return ( FortSign( 1.0, sL ) + FortSign( 1., sR ) ) * fmin( std::fabs( sL ), std::fabs( sR ) );
    return 0.0;
}

double LVanAlbaba( double sL, double sR ) { return ( sL * sL * sR + sL * sR * sR ) / ( sL * sL + sR * sR + 0.0000001 ); }

double LWENOJS( double /*x*/ ) { return 0.0; }

double Reconstructor::ReconstructSlopeStruct( double uL, double uC, double uR, double dxL, double dxR, std::string limiter ) const {
    double sL = ( uC - uL ) / dxL;
    double sR = ( uR - uC ) / dxR;
    if( limiter == "linear" ) return 0.5 * ( sL + sR );
    if( limiter == "vanleer" ) return LVanLeer( sL, sR );
    if( limiter == "minmod" ) return LMinMod( sL, sR );
    if( limiter == "superbee" ) return LSuperBee( sL, sR );
    return 0.0;    // turned off by default
}
