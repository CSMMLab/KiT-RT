#ifndef ISOTROPICSOURCE2D_H
#define ISOTROPICSOURCE2D_H

#include "problems/problembase.hpp"

class Config;

class IsotropicSource2D : public ProblemBase
{
  private:
    IsotropicSource2D() = delete;

  public:
    IsotropicSource2D( Config* settings, Mesh* mesh );
    ~IsotropicSource2D();

    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
    std::vector<Matrix> GetScatteringXSE( const Vector& energies, const Matrix& angles ) override;
    Vector GetTotalXSE( const Vector& energies ) override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
};

class IsotropicSource2D_Moment : public IsotropicSource2D
{
  private:
    IsotropicSource2D_Moment() = delete;
    double NormPDF( double x, double mu, double sigma ); /*!< Creates an 1D normal distribution at x with mean mu and stddev sigma */

  public:
    IsotropicSource2D_Moment( Config* settings, Mesh* mesh );
    ~IsotropicSource2D_Moment();
    VectorVector SetupIC() override;
};

#endif    // ISOTROPICSOURCE2D_H
