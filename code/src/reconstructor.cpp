#include "reconstructor.h"

Reconstructor::Reconstructor(Settings* settings)
{

}

double Reconstructor::Slope(const Vector& Omega, double psiL, double psiR, const Vector& n)const{

    return 0.0;
}

double Sign(double x){

    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;
    return 0.0;
}

double FortSign(double a, double b){

    return abs(a) * Sign(b);
}

double LMinMod(double sL, double sR){

    return 0.5 * (FortSign(1.0, sL) + FortSign(1., sR)) * fmin(abs(sL), abs(sR));

}

double LvanLeer(double sL, double sR){

    return (FortSign(1.0, sL) + FortSign(1.0, sR)) * abs(sL) * abs(sR) / (abs(sL) + abs(sR) + 0.0000001);

}

double LSuperBee(double sL, double sR){

    if (sR >= 0.5 * sL && sR <= 2.0 * sL) return 0.5 * (FortSign(1.0, sL) + FortSign(1., sR)) * fmax(abs(sL), abs(sR));
    if (sR < 0.5 * sL && sR > 2.0 * sL) return (FortSign(1.0, sL) + FortSign(1., sR)) * fmin(abs(sL), abs(sR));
    return 0.0;

}

double LVanAlbaba(double sL, double sR){

    return (sL* sL * sR + sL * sR * sR) / (sL * sL + sR * sR);

}

double LWENOJS(double x){

    return 0.0;

}