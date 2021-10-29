#ifndef WATERPHANTOM_H
#define WATERPHANTOM_H

#include "electronrt.hpp"

class WaterPhantom : public ElectronRT
{
  private:
    WaterPhantom() = delete;

  public:
    WaterPhantom( Config* settings, Mesh* mesh );
    virtual ~WaterPhantom();
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
    std::vector<double> GetDensity( const VectorVector& cellMidPoints );
};

#endif    // WATERPHANTOM_H
