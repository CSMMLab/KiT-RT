#ifndef PHANTOM2D
#define PHANTOM2D

#include "electronrt.hpp"
#include <fstream>
#include <numeric>

#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "toolboxes/interpolation.hpp"

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
