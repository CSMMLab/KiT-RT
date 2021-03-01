#ifndef SLABGEOHG_H
#define SLABGEOHG_H

#include <vector>    // for vector

#include "common/typedef.h"    // for VectorVector, Vector
#include "problembase.h"

class Config;
class Mesh;

class SlabGeoHG : public ProblemBase
{
  private:
    SlabGeoHG() = delete;

  public:
    SlabGeoHG( Config* settings, Mesh* mesh );
    virtual ~SlabGeoHG();

    virtual VectorVector GetScatteringXS( const Vector& energies );
    virtual VectorVector GetTotalXS( const Vector& energies );
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
};

#endif    // SLABGEOHG_H
