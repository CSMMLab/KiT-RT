#ifndef AIRCAVITY
#define AIRCAVITY

#include <vector>    // for vector

#include "common/typedef.h"    // for VectorVector, Vector
#include "electronrt.h"

class Config;
class Mesh;

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
