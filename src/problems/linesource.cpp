#include "problems/linesource.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"
#include "quadratures/quadraturebase.hpp"
#include "velocitybasis/sphericalbase.hpp"
#include "velocitybasis/sphericalharmonics.hpp"
#include <complex>

// ---- Linesource ----

LineSource::LineSource( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _sigmaS = settings->GetSigmaS(); }

LineSource::~LineSource() {}

double LineSource::GetAnalyticalSolution( double x, double y, double t, double /*sigma_s*/ ) {

    double solution = 0.0;
    double R        = sqrt( x * x + y * y );

    if( t > 0 ) {
        if( _sigmaS == 0.0 ) {
            if( ( t - R ) > 0 ) {
                solution = 1 / ( 2 * M_PI * t * sqrt( t * t - R * R ) );
            }
        }
        else if( _sigmaS == 1.0 ) {
            double gamma = R / t;

            if( ( 1 - gamma ) > 0 ) {
                solution = exp( -t ) / ( 2 * M_PI * t * t * sqrt( 1 - gamma * gamma ) ) + 2 * t * HelperIntRho_ptc( R, t );
            }
        }
        else {
            solution = 0.0;
        }
    }

    return M_PI * solution;    // Scaling of soolution ( 4 * M_PI ) *
}

double LineSource::HelperIntRho_ptc( double R, double t ) {

    int numsteps    = 100;
    double integral = 0;
    double gamma    = R / t;
    double omega    = 0;
    // integral is from 0 to  sqrt( 1 - gamma * gamma )
    double stepsize = sqrt( 1 - gamma * gamma ) / (double)numsteps;

    for( int i = 0; i < numsteps; i++ ) {
        omega = i * stepsize + 0.5 * stepsize;

        integral += stepsize * HelperRho_ptc( t * sqrt( gamma * gamma + omega * omega ), t );
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
    double gamma  = R / t;
    double result = 0;
    if( 1 - gamma > 0 ) {
        // Compute the integralpart with midpoint rule
        result = exp( -t ) / ( 32 * M_PI * M_PI * R ) * ( 1 - gamma * gamma ) * HelperIntRho_ptc2( t, gamma );
    }
    return result;
}

double LineSource::HelperIntRho_ptc2( double t, double gamma ) {
    double u        = 0;
    double integral = 0;
    std::complex<double> beta( 0, 0 );
    double q = 0;
    std::complex<double> complexPart( 0, 0 );
    std::complex<double> com_one( 0, 1 );

    // compute the integral using the midpoint rule
    // integral is from 0 to pi

    int numsteps = 100;

    double stepsize = M_PI / (double)numsteps;

    q = ( 1.0 + gamma ) / ( 1.0 - gamma );

    for( int i = 0; i < numsteps; i++ ) {
        u = i * stepsize + 0.5 * stepsize;

        // function evaluation
        beta        = ( log( q ) + com_one * u ) / ( gamma + com_one * tan( 0.5 * u ) );
        complexPart = ( gamma + com_one * tan( 0.5 * u ) ) * beta * beta * beta * exp( 0.5 * t * ( 1 - gamma * gamma ) * beta );
        integral += stepsize * ( 1 / cos( 0.5 * u ) ) * ( 1 / cos( 0.5 * u ) ) * complexPart.real();
    }
    return integral;
}

// ---- LineSource_SN ----

LineSource_SN::LineSource_SN( Config* settings, Mesh* mesh ) : LineSource( settings, mesh ) {}

LineSource_SN::~LineSource_SN() {}

VectorVector LineSource_SN::GetScatteringXS( const Vector& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) );
}

VectorVector LineSource_SN::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) ); }

std::vector<VectorVector> LineSource_SN::GetExternalSource( const Vector& /*energies*/ ) {
    return std::vector<VectorVector>( 1u, std::vector<Vector>( _mesh->GetNumCells(), Vector( 1u, 0.0 ) ) );
}

VectorVector LineSource_SN::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    double t      = 3.2e-4;    // pseudo time for gaussian smoothing
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];
        psi[j]   = Vector( _settings->GetNQuadPoints(), 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) ) ) / ( 4 * M_PI );
    }
    return psi;
}

// ---- LineSource_PN ----

LineSource_Moment::LineSource_Moment( Config* settings, Mesh* mesh ) : LineSource( settings, mesh ) {}

LineSource_Moment::~LineSource_Moment() {}

VectorVector LineSource_Moment::GetScatteringXS( const Vector& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) );
}

VectorVector LineSource_Moment::GetTotalXS( const Vector& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) );
}

std::vector<VectorVector> LineSource_Moment::GetExternalSource( const Vector& /*energies*/ ) {
    SphericalBase* tempBase  = SphericalBase::Create( _settings );
    unsigned ntotalEquations = tempBase->GetBasisSize();
    delete tempBase;    // Only temporally needed

    return std::vector<VectorVector>( 1u, std::vector<Vector>( _mesh->GetNumCells(), Vector( ntotalEquations, 0.0 ) ) );
}

VectorVector LineSource_Moment::SetupIC() {
    // Compute number of equations in the system

    // In case of PN, spherical basis is per default SPHERICAL_HARMONICS
    SphericalBase* tempBase  = SphericalBase::Create( _settings );
    unsigned ntotalEquations = tempBase->GetBasisSize();

    VectorVector initial_sol( _mesh->GetNumCells(), Vector( ntotalEquations, 0 ) );    // zero could lead to problems?
    VectorVector cellMids = _mesh->GetCellMidPoints();

    Vector uIC( ntotalEquations, 0 );

    if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS || _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS_ROTATED ) {
        QuadratureBase* quad          = QuadratureBase::Create( _settings );
        VectorVector quadPointsSphere = quad->GetPointsSphere();
        Vector w                      = quad->GetWeights();

        double my, phi;
        VectorVector moments = VectorVector( quad->GetNq(), Vector( tempBase->GetBasisSize(), 0.0 ) );

        for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
            my                = quadPointsSphere[idx_quad][0];
            phi               = quadPointsSphere[idx_quad][1];
            moments[idx_quad] = tempBase->ComputeSphericalBasis( my, phi );
        }
        // Integrate <1*m> to get factors for monomial basis in isotropic scattering
        for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
            uIC += w[idx_quad] * moments[idx_quad];
        }
        delete quad;
    }

    // Initial condition is dirac impulse at (x,y) = (0,0) ==> constant in angle ==> all moments - exept first - are zero.
    double t       = 3.2e-4;    // pseudo time for gaussian smoothing (Approx to dirac impulse)
    double epsilon = 1e-4;      // minimal value for first moment to avoid div by zero error

    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];    // (x- 0.5) * (x- 0.5)

        double kinetic_density = std::max( 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) ), epsilon );

        if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS || _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS_ROTATED ) {
            initial_sol[j] = kinetic_density * uIC / uIC[0];    // Remember scaling
        }
        if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
            initial_sol[j][0] = kinetic_density;
        }
    }
    delete tempBase;    // Only temporally needed
    return initial_sol;
}

// ---- LineSource SN pseudo1D ----

LineSource_SN_1D::LineSource_SN_1D( Config* settings, Mesh* mesh ) : LineSource_SN( settings, mesh ) {}

VectorVector LineSource_SN_1D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids  = _mesh->GetCellMidPoints();
    double t       = 3.2e-4;    // pseudo time for gaussian smoothing
    double epsilon = 1e-3;      // minimal value for first moment to avoid div by zero error

    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x               = cellMids[j][0];
        double kinetic_density = std::max( 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x ) / ( 4 * t ) ), epsilon );
        psi[j]                 = kinetic_density;
    }
    return psi;
}

// ---- LineSource Moment pseudo1D ----

LineSource_Moment_1D::LineSource_Moment_1D( Config* settings, Mesh* mesh ) : LineSource_Moment( settings, mesh ) {}

VectorVector LineSource_Moment_1D::SetupIC() {
    double t       = 3.2e-4;    // pseudo time for gaussian smoothing (Approx to dirac impulse)
    double epsilon = 1e-3;      // minimal value for first moment to avoid div by zero error

    // In case of PN, spherical basis is per default SPHERICAL_HARMONICS
    if( _settings->GetSolverName() == PN_SOLVER || _settings->GetSolverName() == CSD_PN_SOLVER ) {
        // In case of PN, spherical basis is per default SPHERICAL_HARMONICS in 3 velocity dimensions
        SphericalHarmonics* tempBase = new SphericalHarmonics( _settings->GetMaxMomentDegree(), 3 );
        unsigned ntotalEquations     = tempBase->GetBasisSize();
        delete tempBase;
        VectorVector initialSolution( _mesh->GetNumCells(), Vector( ntotalEquations, 0.0 ) );    // zero could lead to problems?
        VectorVector cellMids = _mesh->GetCellMidPoints();

        for( unsigned idx_cell = 0; idx_cell < cellMids.size(); ++idx_cell ) {
            double x                     = cellMids[idx_cell][0];
            double kinetic_density       = std::max( 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x ) / ( 4 * t ) ), epsilon );
            initialSolution[idx_cell][0] = kinetic_density;
        }
        return initialSolution;
    }
    else {
        SphericalBase* tempBase  = SphericalBase::Create( _settings );
        unsigned ntotalEquations = tempBase->GetBasisSize();

        VectorVector initialSolution( _mesh->GetNumCells(), Vector( ntotalEquations, 0 ) );    // zero could lead to problems?
        VectorVector cellMids = _mesh->GetCellMidPoints();
        Vector uIC( ntotalEquations, 0 );

        if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS || _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS_ROTATED ) {
            QuadratureBase* quad          = QuadratureBase::Create( _settings );
            VectorVector quadPointsSphere = quad->GetPointsSphere();
            Vector w                      = quad->GetWeights();

            double my, phi;
            VectorVector moments = VectorVector( quad->GetNq(), Vector( tempBase->GetBasisSize(), 0.0 ) );

            for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
                my                = quadPointsSphere[idx_quad][0];
                phi               = quadPointsSphere[idx_quad][1];
                moments[idx_quad] = tempBase->ComputeSphericalBasis( my, phi );
            }
            // Integrate <1*m> to get factors for monomial basis in isotropic scattering
            for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
                uIC += w[idx_quad] * moments[idx_quad];
            }
            delete quad;
        }
        // Initial condition is dirac impulse at (x,y) = (0,0) ==> constant in angle ==> all moments - exept first - are zero.
        for( unsigned j = 0; j < cellMids.size(); ++j ) {
            double x               = cellMids[j][0];
            double kinetic_density = std::max( 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x ) / ( 4 * t ) ), epsilon );
            if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS  || _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS_ROTATED) {
                initialSolution[j] = kinetic_density * uIC / uIC[0];    // Remember scaling
            }
            if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
                initialSolution[j][0] = kinetic_density;
            }
        }
        delete tempBase;    // Only temporally needed
        return initialSolution;
    }
}
