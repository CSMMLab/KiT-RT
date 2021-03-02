#ifndef SLABGEOHG_H
#define SLABGEOHG_H

#include "common/pch.h"

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
