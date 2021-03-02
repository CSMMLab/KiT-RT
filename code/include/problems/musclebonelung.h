#ifndef MUSCLEBONELUNG
#define MUSCLEBONELUNG

#include "common/pch.h"
#include "electronrt.h"

class Config;
class Mesh;

class MuscleBoneLung : public ElectronRT
{
  private:
    MuscleBoneLung() = delete;

  public:
    MuscleBoneLung( Config* settings, Mesh* mesh );
    virtual ~MuscleBoneLung();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
    std::vector<double> GetDensity( const VectorVector& cellMidPoints );
};

#endif    // MUSCLEBONELUNG
