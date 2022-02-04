#ifndef AIRCAVITY
#define AIRCAVITY

#include "problems/problembase.hpp"

class AirCavity1D : public ProblemBase
{
  private:
    AirCavity1D() = delete;

  public:
    AirCavity1D( Config* settings, Mesh* mesh );
    ~AirCavity1D();
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    VectorVector SetupIC() override;
    std::vector<double> GetDensity( const VectorVector& cellMidPoints ) override;
    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
    std::vector<Matrix> GetScatteringXSE( const Vector& energies, const Matrix& angles ) override;
};

#endif    // AIRCAVITY
