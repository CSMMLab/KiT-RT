#ifndef CUBIC2DSPLINE_H
#define CUBIC2DSPLINE_H

#include "settings/typedef.h"

class Cubic2DSpline
{
  private:
    Vector _x;
    Vector _y;
    Matrix _data;

    /*
     * calculates a cubic 1D interpolation at point x
     * @param[in]: double param[4] - contains the 4 points needed for the purpose of cubic interpolation
     * @param[in]: double x - the desired interpolation point
     * @param[out]: double - interpolation result
     */
    inline double interpolate1D( double param[4], double x );

    /*
     * finds index of the closest distance to provided value in a given vector
     * @param[in]: double value - value to search for
     * @param[in]: Vector v - sorted vector holding the data
     * @param[out]: unsigned - resulting index of closest value in v
     */
    unsigned indexOfClosestValue( double value, const Vector& v );

    Cubic2DSpline() = delete;

  public:
    /*
     * calculates a cubic 2D interpolation at point (paramX, paramY)
     * @param[in]: double x - x coordinate of the desired interpolation point
     * @param[in]: double y - y coordinate of the desired interpolation point
     * @param[out]: double - interpolation result
     */
    double operator()( double x, double y );

    Cubic2DSpline( const Vector& x, const Vector& y, const Matrix& data );
    ~Cubic2DSpline();
};

#endif
