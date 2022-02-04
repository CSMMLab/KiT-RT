#ifndef WATERPHANTOM_H
#define WATERPHANTOM_H

#include "problems/problembase.hpp"

class WaterPhantom1D : public ProblemBase
{
  private:
    WaterPhantom1D() = delete;

  public:
    WaterPhantom1D( Config* settings, Mesh* mesh );
    virtual ~WaterPhantom1D();
    VectorVector SetupIC() override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
    std::vector<Matrix> GetScatteringXSE( const Vector& energies, const Matrix& angles ) override;
};

class WaterPhantom : public ProblemBase
{
  private:
    WaterPhantom() = delete;

  public:
    WaterPhantom( Config* settings, Mesh* mesh );
    virtual ~WaterPhantom();
    VectorVector SetupIC() override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
    std::vector<Matrix> GetScatteringXSE( const Vector& energies, const Matrix& angles ) override;
};

#endif    // WATERPHANTOM_H
