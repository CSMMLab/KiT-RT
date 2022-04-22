#ifndef ISOTROPICSOURCE2D_CT_H
#define ISOTROPICSOURCE2D_CT_H

#include "problems/problembase.hpp"

class RadiationCTImage : public ProblemBase
{ 
  protected:
  double NormPDF( double x, double mu, double sigma ); /*!< Creates an 1D normal distribution at x with mean mu and stddev sigma */
  
  private:
    RadiationCTImage() = delete;

  public:
    RadiationCTImage( Config* settings, Mesh* mesh );
    virtual ~RadiationCTImage();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
    std::vector<double> GetDensity( const VectorVector& cellMidPoints );
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
};

class RadiationCTImage_Moment : public ProblemBase
{
  private:
    RadiationCTImage_Moment() = delete;

  public:
    RadiationCTImage_Moment( Config* settings, Mesh* mesh );
    virtual ~RadiationCTImage_Moment();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
    std::vector<double> GetDensity( const VectorVector& cellMidPoints );
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
};

#endif    // ISOTROPICSOURCE2D_CT_H
