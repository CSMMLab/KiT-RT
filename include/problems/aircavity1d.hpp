#ifndef AIRCAVITY
#define AIRCAVITY

#include "problems/problembase.hpp"

class AirCavity1D : public ProblemBase
{
  private:
    AirCavity1D() = delete;
    double _sigmaS; /*!< @brief Scattering coefficient */

  public:
    AirCavity1D( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~AirCavity1D();
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    VectorVector SetupIC() override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
};

// Moment solver version
class AirCavity1D_Moment : public ProblemBase
{
  private:
    AirCavity1D_Moment() = delete;
    double _sigmaS; /*!< @brief Scattering coefficient */

  public:
    AirCavity1D_Moment( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~AirCavity1D_Moment();
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    VectorVector SetupIC() override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
};

#endif    // AIRCAVITY
