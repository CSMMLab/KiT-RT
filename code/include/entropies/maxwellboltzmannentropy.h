#ifndef MAXWELLBOLTZMANNENTROPY_H
#define MAXWELLBOLTZMANNENTROPY_H

#include "common/pch.h"

class MaxwellBoltzmannEntropy : public EntropyBase
{
  public:
    inline MaxwellBoltzmannEntropy() : EntropyBase() {}

    inline ~MaxwellBoltzmannEntropy() {}

    inline double Entropy( double z ) override { return z * log( z ) - z; }

    inline double EntropyPrime( double z ) override { return log( z ); }

    // inline void EntropyPrime( Vector& alpha, Vector& m, Vector& grad ) override { grad = m * log( dot( alpha, m ) ); }

    inline double EntropyDual( double y ) override { return exp( y ); }

    // inline void EntropyPrimeDual( Vector& alpha, Vector& m, Vector& grad ) override { grad = m * exp( dot( alpha, m ) ); }

    inline double EntropyPrimeDual( double y ) override { return exp( y ); }

    inline double EntropyHessianDual( double y ) override { return exp( y ); }

    inline bool CheckDomain( double z ) override { return ( z >= 0 && std::isnan( z ) ); }
};

#endif    // MAXWELLBOLTZMANNENTROPY_H
