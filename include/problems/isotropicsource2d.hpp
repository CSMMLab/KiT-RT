#ifndef ISOTROPICSOURCE2D_H
#define ISOTROPICSOURCE2D_H

#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "electronrt.hpp"
#include "toolboxes/interpolation.hpp"

#include <fstream>
#include <numeric>

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

class IsotropicSource2D_Moment : public ElectronRT
{
  private:
    IsotropicSource2D_Moment() = delete;

  public:
    IsotropicSource2D_Moment( Config* settings, Mesh* mesh );
    virtual ~IsotropicSource2D_Moment();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
    std::vector<double> GetDensity( const VectorVector& cellMidPoints );
};

#endif    // ISOTROPICSOURCE2D_H
