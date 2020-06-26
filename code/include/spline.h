#ifndef SPLINE_H
#define SPLINE_H

#include <blaze/math/lapack/posv.h>

#include "settings/typedef.h"
#include "toolboxes/errormessages.h"

class Spline
{
  public:
    enum bd_type { first_deriv = 1, second_deriv = 2 };

  private:
    Vector m_x, m_y;
    Vector m_a, m_b, m_c;
    double m_b0, m_c0;
    bd_type m_left, m_right;
    double m_left_value, m_right_value;
    bool m_force_linear_extrapolation;

  public:
    Spline();
    void set_boundary( bd_type left, double left_value, bd_type right, double right_value, bool force_linear_extrapolation = false );
    void set_points( const std::vector<double>& x, const std::vector<double>& y, bool cubic_Spline = true );
    void set_points( const Vector& x, const Vector& y, bool cubic_Spline = true );
    double operator()( double x ) const;
};

#endif
