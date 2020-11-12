#ifndef PHANTOM2D
#define PHANTOM2D

#include "electronrt.h"
#include <fstream>
#include <numeric>

#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "interpolation.h"

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