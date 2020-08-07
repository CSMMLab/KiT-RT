#ifndef LINESOURCE_H
#define LINESOURCE_H

#include "problembase.h"

class LineSource_SN : public ProblemBase
{
  private:
    LineSource_SN() = delete;

  public:
    LineSource_SN( Config* settings, Mesh* mesh );
    virtual ~LineSource_SN();

    virtual VectorVector GetScatteringXS( const Vector& energies );
    virtual VectorVector GetTotalXS( const Vector& energies );
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies );
    virtual VectorVector SetupIC();
};

class LineSource_SN_Pseudo1D : public LineSource_SN
{
  private:
    LineSource_SN_Pseudo1D() = delete;

  public:
    LineSource_SN_Pseudo1D( Config* settings, Mesh* mesh );

    VectorVector SetupIC() override;
};

class LineSource_SN_Pseudo1D_Physics : public LineSource_SN_Pseudo1D
{
  private:
    LineSource_SN_Pseudo1D_Physics() = delete;

  public:
    LineSource_SN_Pseudo1D_Physics( Config* settings, Mesh* mesh );
};

class LineSource_PN : public ProblemBase
{
  private:
    LineSource_PN() = delete;

    /**
     * @brief Gets the global index for given order l of Legendre polynomials and given
     *        order k of Legendre functions.
     *        Note: This is code doubling from PNSolver::GlobalIndex
     * @param l : order of Legendre polynomial
     * @param k : order of Legendre function
     * @returns global index
     */
    int GlobalIndex( int l, int k ) const;

  public:
    LineSource_PN( Config* settings, Mesh* mesh );
    virtual ~LineSource_PN();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
};

#endif    // LINESOURCE_H
