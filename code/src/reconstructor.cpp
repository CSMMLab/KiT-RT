#include "reconstructor.h"
#include "typedef.h"

Reconstructor::Reconstructor( Config* settings ) {}

double Sign( double x ) {

    if( x > 0.0 ) return 1.0;
    if( x < 0.0 ) return -1.0;
    return 0.0;
}

double FortSign( double a, double b ) { return abs( a ) * Sign( b ); }

double LMinMod( double sL, double sR ) { return 0.5 * ( FortSign( 1.0, sL ) + FortSign( 1., sR ) ) * fmin( abs( sL ), abs( sR ) ); }

double LVanLeer( double sL, double sR ) {

    return ( FortSign( 1.0, sL ) + FortSign( 1.0, sR ) ) * abs( sL ) * abs( sR ) / ( abs( sL ) + abs( sR ) + 0.0000001 );
}

double LSuperBee( double sL, double sR ) {

    if( sR >= 0.5 * sL && sR <= 2.0 * sL ) return 0.5 * ( FortSign( 1.0, sL ) + FortSign( 1., sR ) ) * fmax( abs( sL ), abs( sR ) );
    if( sR < 0.5 * sL && sR > 2.0 * sL ) return ( FortSign( 1.0, sL ) + FortSign( 1., sR ) ) * fmin( abs( sL ), abs( sR ) );
    return 0.0;
}

double LVanAlbaba( double sL, double sR ) { return ( sL * sL * sR + sL * sR * sR ) / ( sL * sL + sR * sR + 0.0000001 ); }

double LWENOJS( double x ) { return 0.0; }

double Reconstructor::ReconstructSlopeStruct( double uL, double uC, double uR, double dxL, double dxR, std::string limiter ) const {
    double sL = ( uC - uL ) / dxL;
    double sR = ( uR - uC ) / dxR;
    if( limiter == "linear" ) return 0.5 * ( sL + sR );
    if( limiter == "vanleer" ) return LVanLeer( sL, sR );
    if( limiter == "minmod" ) return LMinMod( sL, sR );
    if( limiter == "superbee" ) return LSuperBee( sL, sR );
    return 0.0;    // turned off by default
}

void Reconstructor::ReconstructSlopeUnstruct( unsigned nq,
                                              unsigned ncell,
                                              VectorVector& psi,
                                              VectorVector& psiDerX,
                                              VectorVector& psiDerY,
                                              Vector& area,
                                              VectorVector& neighbor,
                                              VectorVector& nx,
                                              VectorVector& ny ) const {

    Vector dPsiMax = Vector( nq, 0.0 );
    Vector dPsiMin = Vector( nq, 0.0 );

    // TODO: boundary treatment
    for( unsigned k = 0; k < nq; ++k ) {
        for( unsigned j = 0; j < ncell; ++j ) {
            double phi;
            VectorVector psiSample = std::vector( nq, Vector( neighbor[j].size(), 0.0 ) );
            VectorVector phiSample = std::vector( nq, Vector( neighbor[j].size(), 0.0 ) );

            for( unsigned l = 0; l < neighbor[j].size(); ++l ) {
                // step 1: calculate deltapsi around neighbors
                if( psi[neighbor[j][l]][k] - psi[j][k] > dPsiMax[k] ) dPsiMax[k] = psi[neighbor[j][l]][k] - psi[j][k] > dPsiMax[k];
                if( psi[neighbor[j][l]][k] - psi[j][k] < dPsiMin[k] ) dPsiMin[k] = psi[neighbor[j][l]][k] - psi[j][k] > dPsiMax[k];

                // step 2: choose sample points (now central points)
                psiSample[k][l] = 0.5 * ( psi[j][k] + psi[neighbor[j][l]][k] );

                // step 3: calculate Phi_ij at sample points
                if( psiSample[k][l] > psi[j][k] ) {
                    phiSample[k][l] = fmin( 1.0, dPsiMax[k] / ( psiSample[k][l] - psi[j][k] ) );
                }
                else if( psiSample[l][k] < psi[j][k] ) {
                    phiSample[k][l] = fmin( 1.0, dPsiMin[k] / ( psiSample[k][l] - psi[j][k] ) );
                }
                else {
                    phiSample[k][l] = 1.0;
                }
            }

            // step 4: find minimum limiter function Phi
            phi = min( phiSample[k] );

            // step 5: reconstruct slopes via Gauss theorem
            for( unsigned l = 0; l < neighbor[j].size(); ++l ) {
                psiDerX[j][k] += phi * 0.5 * ( psi[j][k] + psi[neighbor[j][l]][k] ) * nx[j][l] / area[j];
                psiDerY[j][k] += phi * 0.5 * ( psi[j][k] + psi[neighbor[j][l]][k] ) * ny[j][l] / area[j];
            }
        }
    }
}
