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

    virtual VectorVector GetScatteringXS( const std::vector<double>& energies );
    virtual VectorVector GetTotalXS( const std::vector<double>& energies );
    virtual std::vector<VectorVector> GetExternalSource( const std::vector<double>& energies );
    virtual Vector GetStoppingPower( const Vector& energies );
    virtual VectorVector SetupIC();
};

#endif    // CHECKERBOARD_H
