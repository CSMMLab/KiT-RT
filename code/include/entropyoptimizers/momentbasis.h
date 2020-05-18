#ifndef MOMENTBASIS_H
#define MOMENTBASIS_H

class MomentBasis
{
  public:
    MomentBasis() {}

    static double p_0( double v ) { return 1; }
    static double p_1( double v ) { return v; }
    static double p_2( double v ) { return 0.5 * v * v; }
};

#endif    // MOMENTBASIS_H
