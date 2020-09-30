#ifndef LINESOURCE_H
#define LINESOURCE_H

#include "problembase.h"

class LineSource : public ProblemBase
{
  private:
    LineSource() = delete;

  public:
    LineSource( Config* settings, Mesh* mesh );

    ~LineSource();

    /*!< @brief: Exact analytical solution for the Line Source Test Case at
         @param: x: x coordinate of exact solution
                 y: y coordinate of exact solution
                 t: time of the exact solution
                 sigma_s: scattering cross section of the exact solution
         @return: exact solution at x,y,t,scatteringXS
    */
    double GetAnalyticalSolution( double x, double y, double t, double sigma_s ) override;

  private:
    /*!< @brief: Helper Functions to compute the analytic solution for sigma != 0
     * (See publication: Garret,Hauck; Momentum Closures for Linear Kinetic Transport Equations)
         @param: R = distance to origin
         @param: t = time
         @param: sigma = R/t
     */
    double HelperIntRho_ptc( double R, double t );
    double HelperRho_ptc( double R, double t );
    double HelperRho_ptc1( double R, double t );
    double HelperRho_ptc2( double R, double t );
    double HelperIntRho_ptc2( double t, double gamma );
};

class LineSource_SN : public LineSource
{
  private:
    LineSource_SN() = delete;

  public:
    LineSource_SN( Config* settings, Mesh* mesh );
    ~LineSource_SN();

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

    VectorVector GetScatteringXSE( const Vector& energies, const Vector& angles ) override;
    Vector GetTotalXSE( const Vector& energies ) override;
};

class LineSource_PN : public LineSource
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
    ~LineSource_PN();

    virtual VectorVector GetScatteringXS( const Vector& energies ) override;
    virtual VectorVector GetTotalXS( const Vector& energies ) override;
    virtual std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;
    virtual VectorVector SetupIC() override;
};

#endif    // LINESOURCE_H
