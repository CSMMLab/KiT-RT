#ifndef CHECKERBOARD_H
#define CHECKERBOARD_H

#include "problembase.h"

class Checkerboard : public ProblemBase
{
  private:
    Vector _scatteringXS;
    Vector _totalXS;

    Checkerboard() = delete;

    bool isAbsorption( const Vector& pos );

  public:
    Checkerboard( Config* settings, Mesh* mesh );
    virtual ~Checkerboard();

    virtual VectorVector GetScatteringXS( const std::vector<double>& energies );
    virtual VectorVector GetTotalXS( const std::vector<double>& energies );
    virtual std::vector<double> GetStoppingPower( const std::vector<double>& energies );
    virtual VectorVector SetupIC();
};

#endif    // CHECKERBOARD_H
