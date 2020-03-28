#include "upwindflux.h"

UpwindFlux::UpwindFlux(Settings* settings) : NumericalFlux(settings)
{

}

double UpwindFlux::Flux(const Vector& Omega, double psiL, double psiR, const Vector& n)const{
    /*if( Omega * n > 0 ) {
        return Omega*psiL * norm(n);
    }
    else {
        return Omega*psiR * norm(n);
    }*/
    return -1.0;
}
