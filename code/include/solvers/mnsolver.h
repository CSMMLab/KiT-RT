#ifndef MNSOLVER_H
#define MNSOLVER_H

#include <pnsolver.h>

class MNSolver : PNSolver
{
    /*! @brief: compute the legendre function of degree l and order k  at point x using the lth degree
     *          legendre polynomial
     * @param: [in] double x : spatial point, -1 <= x <= 1
     *         [in] int    l : degree of associated legendre polynomial, 0 <= l
     *         [in] int    k : order of associated legendre polynomial, -l <= k <= l
     *         [ou] double   : value of legendre polynomial of degree l and order k  at point x
     */
    double AssociatedLegendrePoly( double x, int l, int k );
#endif    // MNSOLVER_H
