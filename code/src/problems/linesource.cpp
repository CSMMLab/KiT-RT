#include "problems/linesource.h"
#include "common/config.h"
#include "common/mesh.h"
#include "physics.h"

#include <complex>

// ---- Linesource ----

LineSource::LineSource( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _physics = nullptr; }
LineSource::~LineSource() {}

double LineSource::GetAnalyticalSolution( double x, double y, double t, double sigma_s ) {

    double solution = 0.0;
    double R        = sqrt( x * x + y * y );

    if( sigma_s == 0.0 ) {
        if( ( t - R ) > 0 ) {
            solution = 1 / ( 2 * M_PI * t * sqrt( t * t - R * R ) );
        }
    }
    else {
        double gamma = R / t;

        if( ( 1 - gamma ) > 0 ) {
            solution = exp( -t ) / ( 2 * M_PI * t * t * sqrt( 1 - gamma * gamma ) ) + 2 * t * HelperIntRho_ptc( R, t );
        }
    }
    return solution;
}

double LineSource::HelperIntRho_ptc( double R, double t ) {

    int numsteps    = 100;
    double integral = 0;
    double gamma    = R / t;
    double omega    = 0;
    // integral is from 0 to  sqrt( 1 - gamma * gamma )
    double stepsize = sqrt( 1 - gamma * gamma ) / (double)numsteps;

    omega = 0.5 * stepsize;
    for( int i = 0; i < numsteps; i++ ) {
        omega = omega + i * stepsize;
        integral += HelperRho_ptc( t * sqrt( gamma * gamma + omega * omega ), t );
    }
    return integral;
}

double LineSource::HelperRho_ptc( double R, double t ) {
    double result = HelperRho_ptc1( R, t ) + HelperRho_ptc2( R, t );

    return result;
}

double LineSource::HelperRho_ptc1( double R, double t ) {
    double gamma  = R / t;
    double result = exp( -t ) / ( 4.0 * M_PI * R * t ) * log( ( 1 + gamma ) / ( 1 - gamma ) );

    return result;
}

double LineSource::HelperRho_ptc2( double R, double t ) {
    double gamma = R / t;

    double result   = 0;
    double integral = 0;
    if( 1 - gamma > 0 ) {

        // Compute the integralpart with midpoint rule
        result = exp( -t ) / ( 32 * M_PI * M_PI * R ) * ( 1 - gamma * gamma ) * HelperIntRho_ptc2( t, gamma );
    }

    return result;
}

double LineSource::HelperIntRho_ptc2( double t, double gamma ) {

    double u        = 0;
    double integral = 0;
    double beta     = 0;
    double q        = 0;
    complex<double> complexPart;

    // compute the integral using the midpoint rule
    // integral is from 0 to pi

    int numsteps = 100;

    double stepsize = M_PI / (double)numsteps;

    u = 0.5 * stepsize;
    q = ( 1 + gamma ) / ( 1 - gamma );

    for( int i = 0; i < numsteps; i++ ) {
        u = u + i * stepsize;

        // function evaluation
        beta        = ( log( q ) + _Complex_I * u ) / ( gamma + _Complex_I * tan( 0.5 * u ) );
        complexPart = ( gamma + _Complex_I * tan( 0.5 * u ) ) * beta * beta * beta * exp( 0.5 * t * ( 1 - gamma * gamma ) * beta );
        integral += ( 1 / cos( 0.5 * u ) ) * ( 1 / cos( 0.5 * u ) ) * complexPart.real();
    }

    return integral;
}

// ---- LineSource_SN ----

LineSource_SN::LineSource_SN( Config* settings, Mesh* mesh ) : LineSource( settings, mesh ) {}

LineSource_SN::~LineSource_SN() {}

VectorVector LineSource_SN::GetScatteringXS( const Vector& energies ) { return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 1.0 ) ); }

VectorVector LineSource_SN::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 1.0 ) ); }

std::vector<VectorVector> LineSource_SN::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( 1u, std::vector<Vector>( _mesh->GetNumCells(), Vector( 1u, 0.0 ) ) );
}

VectorVector LineSource_SN::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    double t      = 3.2e-4;    // pseudo time for gaussian smoothing
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];
        psi[j]   = 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) );
    }
    return psi;
}

// ---- LineSource_SN_peudo1D ----

LineSource_SN_Pseudo1D::LineSource_SN_Pseudo1D( Config* settings, Mesh* mesh ) : LineSource_SN( settings, mesh ) {}

VectorVector LineSource_SN_Pseudo1D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    double t      = 3.2e-4;    // pseudo time for gaussian smoothing
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        psi[j]   = 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x ) / ( 4 * t ) );
    }
    return psi;
}

// ---- LineSource_SN_Pseudo1D_Physics ----

LineSource_SN_Pseudo1D_Physics::LineSource_SN_Pseudo1D_Physics( Config* settings, Mesh* mesh ) : LineSource_SN_Pseudo1D( settings, mesh ) {
    _physics = new Physics( settings->GetHydrogenFile(), settings->GetOxygenFile(), "../input/stopping_power.txt" );
}

VectorVector LineSource_SN_Pseudo1D_Physics::GetScatteringXSE( const Vector& energies, const Vector& angles ) {
    return _physics->GetScatteringXS( energies, angles );
}

Vector LineSource_SN_Pseudo1D_Physics::GetTotalXSE( const Vector& energies ) { return _physics->GetTotalXSE( energies ); }

// ---- LineSource_PN ----

int LineSource_PN::GlobalIndex( int l, int k ) const {
    int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
    int prevIndicesThisLevel = k + l;    // number of previous indices in current level
    return numIndicesPrevLevel + prevIndicesThisLevel;
}

LineSource_PN::LineSource_PN( Config* settings, Mesh* mesh ) : LineSource( settings, mesh ) {}

LineSource_PN::~LineSource_PN() {}

VectorVector LineSource_PN::GetScatteringXS( const Vector& energies ) { return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 1.0 ) ); }

VectorVector LineSource_PN::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 1.0 ) ); }

std::vector<VectorVector> LineSource_PN::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( 1u, std::vector<Vector>( _mesh->GetNumCells(), Vector( 1u, 0.0 ) ) );
}

VectorVector LineSource_PN::SetupIC() {
    // Compute number of equations in the system
    int ntotalEquations = GlobalIndex( _settings->GetMaxMomentDegree(), _settings->GetMaxMomentDegree() ) + 1;

    VectorVector psi( _mesh->GetNumCells(), Vector( ntotalEquations, 0 ) );    // zero could lead to problems?
    VectorVector cellMids = _mesh->GetCellMidPoints();

    // Initial condition is dirac impulse at (x,y) = (0,0) ==> constant in angle ==> all moments are zero.
    double t = 3.2e-4;    // pseudo time for gaussian smoothing (Approx to dirac impulse)
    // double offset = 0.099;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];    // (x- 0.5) * (x- 0.5)

        psi[j][0] = 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( ( x ) * ( x ) + ( y ) * ( y ) ) / ( 4 * t ) );
        //+
        //  1.0 / ( 4.0 * M_PI * t ) * std::exp( -( ( x - offset ) * ( x - offset ) + ( y - offset ) * ( y - offset ) ) / ( 4 * t
        //  ) ) + 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( ( x + offset ) * ( x + offset ) + ( y - offset ) * ( y - offset ) ) / (
        //  4 * t ) ) + 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( ( x - offset ) * ( x - offset ) + ( y + offset ) * ( y + offset )
        //  ) / ( 4 * t ) ) + 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( ( x + offset ) * ( x + offset ) + ( y + offset ) * ( y +
        //  offset ) ) / ( 4 * t ) ) + 1.0 / ( 4.0 * M_PI * t ) *
        //      std::exp( -( ( x - 2 * offset ) * ( x - 2 * offset ) + ( y - 2 * offset ) * ( y - 2 * offset ) ) / ( 4 * t ) ) +
        //  1.0 / ( 4.0 * M_PI * t ) *
        //      std::exp( -( ( x + 2 * offset ) * ( x + 2 * offset ) + ( y - 2 * offset ) * ( y - 2 * offset ) ) / ( 4 * t ) ) +
        //  1.0 / ( 4.0 * M_PI * t ) *
        //      std::exp( -( ( x - 2 * offset ) * ( x - 2 * offset ) + ( y + 2 * offset ) * ( y + 2 * offset ) ) / ( 4 * t ) ) +
        //  1.0 / ( 4.0 * M_PI * t ) * std::exp( -( ( x + 2 * offset ) * ( x + 2 * offset ) + ( y + 2 * offset ) * ( y + 2 *
        //  offset ) ) / ( 4 * t ) );
    }

    // Debugging jump test case
    // for( unsigned j = 0; j < cellMids.size(); ++j ) {
    //    if( cellMids[j][0] > -0.2 && cellMids[j][1] > -0.2 && cellMids[j][1] < 0.2 && cellMids[j][0] < 0.2 )
    //        psi[j][0] = 1.0;
    //    else
    //        psi[j][0] = 0.0;
    //}
    return psi;
}
