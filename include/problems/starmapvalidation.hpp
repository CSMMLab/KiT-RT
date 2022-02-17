#ifndef ISOTROPICSOURCE2D_H
#define ISOTROPICSOURCE2D_H

#include "problems/problembase.hpp"

class Config;

class StarMapValidation_SN : public ProblemBase
{
  private:
    StarMapValidation_SN() = delete;

  public:
    StarMapValidation_SN( Config* settings, Mesh* mesh );
    ~StarMapValidation_SN();

    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
    std::vector<Matrix> GetScatteringXSE( const Vector& energies, const Matrix& angles ) override;
    Vector GetTotalXSE( const Vector& energies ) override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
};

class StarMapValidation_Moment : public StarMapValidation_SN
{
  private:
    StarMapValidation_Moment() = delete;
    double NormPDF( double x, double mu, double sigma ); /*!< Creates an 1D normal distribution at x with mean mu and stddev sigma */

  public:
    StarMapValidation_Moment( Config* settings, Mesh* mesh );
    ~StarMapValidation_Moment();
    VectorVector SetupIC() override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
};

#endif    // ISOTROPICSOURCE2D_H
