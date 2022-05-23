#ifndef RADIATION_CT_H
#define RADIATION_CT_H

#include "problems/problembase.hpp"

class RadiationCTImage : public ProblemBase
{
  private:
    RadiationCTImage() = delete;

  public:
    RadiationCTImage( Config* settings, Mesh* mesh );
    ~RadiationCTImage();
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    VectorVector SetupIC() override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override final;
    VectorVector GetScatteringXS( const Vector& energies ) override final;
    VectorVector GetTotalXS( const Vector& energies ) override final;
};

class RadiationCTImage_Moment : public RadiationCTImage
{
  private:
    RadiationCTImage_Moment() = delete;

  public:
    RadiationCTImage_Moment( Config* settings, Mesh* mesh );
    ~RadiationCTImage_Moment();
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override final;
    VectorVector SetupIC() override final;
};

#endif    // RADIATION_CT_H
