#ifndef PHANTOM2D
#define PHANTOM2D

#include "common/pch.h"
#include "electronrt.h"

class Config;
class Mesh;

class Phantom2D : public ElectronRT
{
  private:
    Phantom2D() = delete;

  public:
    Phantom2D( Config* settings, Mesh* mesh );
    virtual ~Phantom2D();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
    std::vector<double> GetDensity( const VectorVector& cellMidPoints );
};

#endif    // PHANTOM2D
