/*!
 * @file interpolation.h
 * @brief Class for interpolation of tabulated values (used for database interpolation)
 */
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
     * @param[in] x - table values for x
     * @param[in] y - table values for y
     * @param[in] type - linear or cubic interpolation
     */
    Interpolation( const std::vector<double>& x, const std::vector<double>& y, TYPE type = linear );

    /*!
     * @brief constructor linear interpolation for Vector input
     * @param[in] x - table values for x
     * @param[in] y - table values for y
     * @param[in] type - linear or cubic interpolation
     */
    Interpolation( const Vector& x, const Vector& y, TYPE type = linear );

    /*!
     * @brief constructor cubic interpolation
     * @param[in] x - table values for x
     * @param[in] y - table values for y
     * @param[in] type - linear or cubic interpolation
     */
    Interpolation( const Vector& x, const Vector& y, const Matrix& data, TYPE type = cubic );

    // Here the interpolation between the previously defined table values is performed
    /*!
     * @brief defines one dimensional interpolation at x
     * @param[in] x - value at which to interpolate
     * @param[out] y - corresponding interpolated y value
     */
    double operator()( double x ) const;

    /*!
     * @brief defines interpolation for a Vector of values
     * @param[in] v - values at which to interpolate
     * @param[out] y - corresponding interpolated values
     */
    Vector operator()( Vector v ) const;

    /*!
     * @brief defines interpolation for a std::vector of values
     * @param[in] v - values at which to interpolate
     * @param[out] y - corresponding interpolated values
     */
    std::vector<double> operator()( std::vector<double> v ) const;

    /*!
     * @brief defines 2D interpolation at x and y
     * @param[in] x - value at which to interpolate
     * @param[in] y - value at which to interpolate
     * @param[out] data - corresponding interpolated value
     */
    double operator()( double x, double y ) const;
};

#endif
