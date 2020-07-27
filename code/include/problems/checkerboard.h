#ifndef CHECKERBOARD_H
#define CHECKERBOARD_H

#include "problembase.h"

class Checkerboard : public ProblemBase
{
  private:
    Vector _scatteringXS;
    Vector _totalXS;

    Checkerboard() = delete;

    bool isAbsorption( const Vector& pos ) const;
    bool isSource( const Vector& pos ) const;

  public:
    Checkerboard( Config* settings, Mesh* mesh );
    virtual ~Checkerboard();

    virtual VectorVector GetScatteringXS( const Vector& energies );
    virtual VectorVector GetTotalXS( const Vector& energies );
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual Vector GetStoppingPower( const Vector& energies );
    virtual VectorVector SetupIC();
};

#endif    // CHECKERBOARD_H
