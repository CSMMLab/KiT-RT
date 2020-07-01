#ifndef MAXWELLBOLTZMANNENTROPY_H
#define MAXWELLBOLTZMANNENTROPY_H

#include "entropybase.h"
#include <cmath>

class MaxwellBoltzmannEntropy : public EntropyBase
{
  public:
    inline MaxwellBoltzmannEntropy() : EntropyBase() {}

    inline ~MaxwellBoltzmannEntropy() {}

    inline double Entropy( double z ) override { return z * log( z ) - z; }

    inline double EntropyPrime( double z ) override { return log( z ); }

    inline double EntropyDual( double y ) override { return exp( y ); }

    inline double EntropyPrimeDual( double y ) override { return exp( y ); }

    inline bool CheckDomain( double z ) override { return ( z >= 0 && std::isnan( z ) ); }
};

#endif    // MAXWELLBOLTZMANNENTROPY_H
