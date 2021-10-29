#ifndef MUSCLEBONELUNG
#define MUSCLEBONELUNG

#include "electronrt.hpp"

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
