#include "reconstructor.h"

Reconstructor::Reconstructor(Settings* settings)
{

}

double Reconstructor::Slope(const Vector& Omega, double psiL, double psiR, const Vector& n)const{

    return 0.0;
}

int FortSign(double x){

    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

double LMinMod(double x){

    return -1.0;

}

double LWarmBeam(double x){

    return -1.0;

}

double LMuscl(double x){

    return -1.0;

}

double LSuperBee(double x){

    return -1.0;

}

double LOsherC(double x){

    return -1.0;

}

double LWENOJS(double x){

    return -1.0;

}