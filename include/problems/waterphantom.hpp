#ifndef WATERPHANTOM_H
#define WATERPHANTOM_H

#include "problems/problembase.hpp"

class WaterPhantom1D : public ProblemBase
{
  private:
    WaterPhantom1D() = delete;
    double _sigmaS; /*!< @brief Scattering coefficient */
  public:
    WaterPhantom1D( Config* settings, Mesh* mesh );
    virtual ~WaterPhantom1D();
    VectorVector SetupIC() override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
};

class WaterPhantom : public ProblemBase
{
  private:
    WaterPhantom() = delete;
    double _sigmaS; /*!< @brief Scattering coefficient */
  public:
    WaterPhantom( Config* settings, Mesh* mesh );
    virtual ~WaterPhantom();
    VectorVector SetupIC() override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
};

#endif    // WATERPHANTOM_H
