#ifndef QUADRATICENTROPY_H
#define QUADRATICENTROPY_H

#include "entropybase.hpp"
#include <cmath>

class QuadraticEntropy : public EntropyBase
{
  public:
    inline QuadraticEntropy() : EntropyBase() {}

    inline ~QuadraticEntropy() {}

    inline double Entropy( double z ) override { return 0.5 * z * z; }

    inline double EntropyPrime( double z ) override { return z; }

    inline double EntropyDual( double y ) override { return 0.5 * y * y; }

    // inline void EntropyPrimeDual( Vector& alpha, Vector& m, Vector& grad ) override { grad = m * dot( alpha, m ); }

    inline double EntropyPrimeDual( double y ) override { return y; }

    inline double EntropyHessianDual( double /*y*/ ) override { return 1.0; }

    inline bool CheckDomain( double z ) override { return std::isnan( z ); }
};

#endif    // QUADRATICENTROPY_H
