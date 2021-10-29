#ifndef AIRCAVITY
#define AIRCAVITY

#include "electronrt.hpp"

class AirCavity1D : public ElectronRT
{
  private:
    AirCavity1D() = delete;

  public:
    AirCavity1D( Config* settings, Mesh* mesh );
    virtual ~AirCavity1D();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
    std::vector<double> GetDensity( const VectorVector& cellMidPoints );
};

#endif    // AIRCAVITY
