/*!
 * @file interpolation.h
 * @brief Class for interpolation of tabulated values (used for database interpolation)
 */
#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "common/typedef.hpp"

class Interpolation
{
  public:
    enum TYPE { linear, loglinear, cubic };

  private:
    unsigned _dim;
    Vector _x;    /*!< @brief x input data */
    Vector _y;    /*!< @brief y input data */
    Matrix _data; /*!< @brief data matrix w.r.t. x and y grid */
    TYPE _type;   /*!< @brief type of interpolation (linear, loglinear, cubic) */

    Interpolation() = delete;

    /*!
     * @brief Check if vector lengths match and values are sorted ascendingly
     */
    void Setup();
    /*!
     * @return index of element closest to 'value' in 'v'
     */
    unsigned IndexOfClosestValue( double value, const Vector& v ) const;
    /*!
     * @return value of third degree/cubic polynomial with parameters 'param' at 'x'
     */
    inline double EvalCubic1DSpline( double param[4], double x ) const;

  public:
    // In constructor tabulated values are initialised
    /*!
     * @brief constructor linear interpolation for std::vector input
     * @param x table values for x
     * @param y table values for y
     * @param type linear or cubic interpolation
     */
    Interpolation( const std::vector<double>& x, const std::vector<double>& y, TYPE type = linear );

    /*!
     * @brief constructor linear interpolation for Vector input
     * @param x table values for x
     * @param y table values for y
     * @param type linear or cubic interpolation
     */
    Interpolation( const Vector& x, const Vector& y, TYPE type = linear );

    /*!
     * @brief constructor cubic interpolation
     * @param x table values for x
     * @param y table values for y
     * @param data matrix w.r.t x y grid
     * @param type of interpolation  (linear, loglinear, cubic)
     */
    Interpolation( const Vector& x, const Vector& y, const Matrix& data, TYPE type = cubic );

    // Here the interpolation between the previously defined table values is performed
    /*!
     * @brief defines one dimensional interpolation at x
     * @param x value at which to interpolate
     * @returns corresponding interpolated y value
     */
    double operator()( double x ) const;

    /*!
     * @brief defines interpolation for a Vector of values
     * @param v values at which to interpolate
     * @returns corresponding interpolated values
     */
    Vector operator()( Vector v ) const;

    /*!
     * @brief defines interpolation for a std::vector of values
     * @param v values at which to interpolate
     * @returns corresponding interpolated values
     */
    std::vector<double> operator()( std::vector<double> v ) const;

    /*!
     * @brief defines 2D interpolation at x and y
     * @param x value at which to interpolate
     * @param y value at which to interpolate
     * @returns corresponding interpolated value
     */
    double operator()( double x, double y ) const;
};

#endif
