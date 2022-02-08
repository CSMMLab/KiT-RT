#ifndef ISOTROPICSOURCE2D_CT_H
#define ISOTROPICSOURCE2D_CT_H

#include "problems/problembase.hpp"

class RadiationCTImage : public ProblemBase
{
  private:
    RadiationCTImage() = delete;

  public:
    RadiationCTImage( Config* settings, Mesh* mesh );
    virtual ~RadiationCTImage();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
    std::vector<double> GetDensity( const VectorVector& cellMidPoints );
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
};

#endif    // ISOTROPICSOURCE2D_CT_H
