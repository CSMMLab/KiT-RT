#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "common/typedef.h"

class Interpolation
{
  public:
    enum BOUNDARY { first_deriv, second_deriv };
    enum TYPE { linear, loglinear, cubic };

  private:
    Vector _x, _y;
    Vector _a, _b, _c;
    double _b0, _c0;
    BOUNDARY _left, _right;
    TYPE _type;
    double _left_value, _right_value;
    bool _force_linear_extrapolation;

    Interpolation() = delete;

  public:
    Interpolation( const std::vector<double>& x, const std::vector<double>& y, TYPE type = linear );
    Interpolation( const Vector& x, const Vector& y, TYPE type = linear );

    void set_boundary( BOUNDARY left, double left_value, BOUNDARY right, double right_value, bool force_linear_extrapolation = false );
    void set_points( const std::vector<double>& x, const std::vector<double>& y, TYPE type );
    void set_points( const Vector& x, const Vector& y, TYPE type );
    double operator()( double x ) const;
    Vector operator()( Vector v ) const;
};

#endif
