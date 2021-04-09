#ifndef ISOTROPICSOURCE2D_H
#define ISOTROPICSOURCE2D_H

#include "electronrt.h"
#include <fstream>
#include <numeric>

#include "common/config.h"
#include "common/io.h"
#include "common/mesh.h"
#include "toolboxes/interpolation.h"

class IsotropicSource2D : public ElectronRT
{
  private:
    IsotropicSource2D() = delete;

  public:
    IsotropicSource2D( Config* settings, Mesh* mesh );
    virtual ~IsotropicSource2D();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
    std::vector<double> GetDensity( const VectorVector& cellMidPoints );
};

#endif    // ISOTROPICSOURCE2D_H
