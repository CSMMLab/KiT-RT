#ifndef CHECKERBOARD_H
#define CHECKERBOARD_H

#include "problembase.h"

class Checkerboard : public ProblemBase
{
  private:
    std::vector<Matrix> _scatteringXS;
    std::vector<double> _totalXS;

    Checkerboard() = delete;

    bool isAbsorption( const Vector& pos );

  public:
    Checkerboard( Config* settings, Mesh* mesh );
    virtual ~Checkerboard();

    virtual std::vector<Matrix> GetScatteringXS( const double energy );
    virtual std::vector<double> GetTotalXS( const double energy );
    virtual std::vector<double> GetStoppingPower( const std::vector<double>& energies );
    virtual VectorVector SetupIC();
};

#endif    // CHECKERBOARD_H
