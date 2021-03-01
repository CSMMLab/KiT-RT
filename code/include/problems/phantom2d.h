#ifndef PHANTOM2D
#define PHANTOM2D

#include "electronrt.h"
#include <vector>    // for vector

#include "common/typedef.h"    // for VectorVector, Vector

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
