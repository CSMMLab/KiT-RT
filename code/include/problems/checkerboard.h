#ifndef CHECKERBOARD_H
#define CHECKERBOARD_H

#include "problembase.h"

class Checkerboard_SN : public ProblemBase
{
  private:
    Vector _scatteringXS;
    Vector _totalXS;

    Checkerboard_SN() = delete;

    bool isAbsorption( const Vector& pos ) const;
    bool isSource( const Vector& pos ) const;

  public:
    Checkerboard_SN( Config* settings, Mesh* mesh );
    virtual ~Checkerboard_SN();

    virtual VectorVector GetScatteringXS( const Vector& energies );
    virtual VectorVector GetTotalXS( const Vector& energies );
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
};

class Checkerboard_PN : public ProblemBase
{
  private:
    Vector _scatteringXS;
    Vector _totalXS;

    Checkerboard_PN() = delete;

    bool isAbsorption( const Vector& pos ) const;
    bool isSource( const Vector& pos ) const;
    int GlobalIndex( int l, int k ) const;

  public:
    Checkerboard_PN( Config* settings, Mesh* mesh );
    virtual ~Checkerboard_PN();

    virtual VectorVector GetScatteringXS( const Vector& energies );
    virtual VectorVector GetTotalXS( const Vector& energies );
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
};

#endif    // CHECKERBOARD_H
