#ifndef LINESOURCE_H
#define LINESOURCE_H

#include "problembase.h"

class LineSource : public ProblemBase
{
  private:
    LineSource() = delete;

  public:
    LineSource( Config* settings, Mesh* mesh );
    virtual ~LineSource();

    virtual VectorVector GetScatteringXS( const std::vector<double>& energies );
    virtual VectorVector GetTotalXS( const std::vector<double>& energies );
    virtual std::vector<VectorVector> GetExternalSource( const std::vector<double>& energies );
    virtual std::vector<double> GetStoppingPower( const std::vector<double>& energies );
    virtual VectorVector SetupIC();
};

#endif    // LINESOURCE_H
