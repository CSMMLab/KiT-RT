#ifndef PNSOLVER_H
#define PNSOLVER_H

#include <cmath>

#include "solver.h"
#include "typedef.h"

#include "settings/config.h"

class PNSolver : public Solver
{
  public:
    /**
     * @brief PNSolver constructor
     * @param settings stores all needed information
     */
    PNSolver( Config* settings );
    /**
     * @brief Solve functions runs main time loop
     */
    virtual void Solve();
    /**
     * @brief Output solution to VTK file
     */
    virtual void Save() const;
    void SolveMPI();    // can be deleted later
};

#endif    // PNSOLVER_H

// class Problem :
//{
//  protected:
//    Settings* _settings;
//    unsigned _nStates;
//    std::shared_ptr<spdlog::logger> _log;
//    Mesh* _mesh;
//
//    // variance vector
//    Vector _sigma;
//
//    Problem() {}
//
//  public:
//    Problem( Settings* settings );
//    static Problem* Create( Settings* settings );
//    virtual ~Problem();
//    virtual void Solve() {}
//    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) = 0;
//    virtual double ComputeDt( const Matrix& u, double dx, unsigned level ) const;
//    virtual Vector IC( const Vector& x, const Vector& xi )     = 0;
//    virtual Vector LoadIC( const Vector& x, const Vector& xi ) = 0;
//    virtual Matrix ExactSolution( double t, const Matrix& x, const Vector& xi ) const;
//    virtual Matrix Source( const Matrix& uQ ) const;
//    virtual Matrix Source( const Matrix& uQ, const Vector& x, double t, unsigned level ) const;
//    virtual void SourceImplicit( Matrix& uQNew, const Matrix& uQTilde, const Matrix& uQ, const Vector& x, double t, unsigned level ) const;
//    virtual Matrix BoundaryFlux( const Matrix& u, const Vector& nUnit, const Vector& n, unsigned level ) const;
//    virtual void SetMesh( Mesh* mesh ) { _mesh = mesh; }
//};

// class PNEquations : public Solver
//{
//  protected:
//    // moment orders for P_N
//    int _N;
//    double _sigmaA;    // absorption coefficient
//    double _sigmaS;    // scattering coefficient
//    double _sigmaT;    // total crossection
//
//    // System Matrix for x, y and z flux
//    Matrix _Ax;
//    Matrix _Ay;
//    Matrix _Az;
//
//    // Roe Matrices
//    Matrix _AbsAx;
//    Matrix _AbsAz;
//
//    // parameter functions for setting up system matrix
//    double AParam( int l, int k ) const;
//    double BParam( int l, int k ) const;
//    double CParam( int l, int k ) const;
//    double DParam( int l, int k ) const;
//    double EParam( int l, int k ) const;
//    double FParam( int l, int k ) const;
//
//    double CTilde( int l, int k ) const;
//    double DTilde( int l, int k ) const;
//    double ETilde( int l, int k ) const;
//    double FTilde( int l, int k ) const;
//
//    // mathematical + index functions
//    int Sgn( int k ) const;
//    int kPlus( int k ) const;
//    int kMinus( int k ) const;
//    virtual int GlobalIndex( int l, int k ) const;
//
//    // function for setting up system matrices
//    void SetupSystemMatrices();
//
//  public:
//    PNEquations( Settings* settings );
//    /**
//     * @brief PNEquations constructur without setting up system matrices used for radiation hydrodynamics
//     * @param settings Settings pointer
//     * @param noSystemMatrix dummy bool
//     */
//    PNEquations( Settings* settings, bool noSystemMatrix );
//    virtual ~PNEquations();
//    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
//    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level );
//    Matrix F( const Vector& u );
//    Matrix F( const Matrix& u );
//    virtual Matrix Source( const Matrix& uQ, const Vector& x, double t, unsigned level ) const;
//    virtual double ComputeDt( const Matrix& u, double dx, unsigned level ) const;
//    virtual Vector IC( const Vector& x, const Vector& xi );
//    virtual Vector LoadIC( const Vector& x, const Vector& xi );
//};

#endif    // PNSOLVER_H
