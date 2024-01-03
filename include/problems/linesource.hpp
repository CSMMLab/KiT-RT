#ifndef LINESOURCE_H
#define LINESOURCE_H

#include "problembase.hpp"

class LineSource : public ProblemBase
{
  private:
    LineSource() = delete;

  protected:
    Vector _sigmaS; /*!< @brief Vector of scattering crosssections */
    Vector _sigmaT; /*!< @brief Vector of total crosssections */

  public:
    LineSource( Config* settings, Mesh* mesh, QuadratureBase* quad );

    ~LineSource();

    /*! @brief Exact analytical solution for the Line Source Test Case at
         @param x coordinate of exact solution
         @param y coordinate of exact solution
         @param t time of the exact solution
         @param sigma_s  scattering cross section of the exact solution
         @return exact solution at x,y,t,scatteringXS
    */
    virtual double GetAnalyticalSolution( double x, double y, double t, double sigma_s ) override;

    VectorVector GetScatteringXS( const Vector& energies ) override;
    VectorVector GetTotalXS( const Vector& energies ) override;
    std::vector<VectorVector> GetExternalSource( const Vector& energies ) override;

  private:
    /*! @brief Helper Functions to compute the analytic solution for sigma != 0
     * (See publication: Garret,Hauck; Momentum Closures for Linear Kinetic Transport Equations)
         @param R distance to origin
         @param t time
     */
    double HelperIntRho_ptc( double R, double t );
    double HelperRho_ptc( double R,
                          double t ); /*!< @brief helper to comput analytic line source solution. @param R distance to origin @param t time */
    double HelperRho_ptc1( double R,
                           double t ); /*!< @brief helper to comput analytic line source solution. @param R distance to origin @param t time */
    double HelperRho_ptc2( double R,
                           double t ); /*!< @brief helper to comput analytic line source solution. @param R distance to origin @param t time */
    double HelperIntRho_ptc2( double t,
                              double gamma ); /*!< @brief helper to comput analytic line source solution @param t time @param gamma equals R/t */
};

class LineSource_SN : public LineSource
{
  private:
    LineSource_SN() = delete;

  public:
    LineSource_SN( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~LineSource_SN();

    VectorVector SetupIC() override;
};

class LineSource_Moment : public LineSource
{
  private:
    LineSource_Moment() = delete;

  public:
    LineSource_Moment( Config* settings, Mesh* mesh, QuadratureBase* quad );
    ~LineSource_Moment();

    VectorVector SetupIC() override;
};

class LineSource_SN_1D : public LineSource_SN
{
  private:
    LineSource_SN_1D() = delete;

  public:
    LineSource_SN_1D( Config* settings, Mesh* mesh, QuadratureBase* quad );

    VectorVector SetupIC() override;
};

class LineSource_Moment_1D : public LineSource_Moment
{
  private:
    LineSource_Moment_1D() = delete;

  public:
    LineSource_Moment_1D( Config* settings, Mesh* mesh, QuadratureBase* quad );

    VectorVector SetupIC() override;
};

#endif    // LINESOURCE_H
