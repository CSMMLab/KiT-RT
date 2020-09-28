#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "common/typedef.h"

class Interpolation
{
  public:
    enum TYPE { linear, loglinear, cubic };

  private:
    unsigned _dim;
    Vector _x, _y;
    Matrix _data;
    TYPE _type;

    Interpolation() = delete;

    void Setup();
    unsigned IndexOfClosestValue( double value, const Vector& v ) const;
    inline double EvalCubic1DSpline( double param[4], double x ) const;

  public:
    Interpolation( const std::vector<double>& x, const std::vector<double>& y, TYPE type = linear );
    Interpolation( const Vector& x, const Vector& y, TYPE type = linear );
    Interpolation( const Vector& x, const Vector& y, const Matrix& data, TYPE type = cubic );

    double operator()( double x ) const;
    Vector operator()( Vector v ) const;
    std::vector<double> operator()( std::vector<double> v ) const;
    double operator()( double x, double y ) const;
};

#endif
