#include "problems/meltingcube.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"
#include "quadratures/quadraturebase.hpp"
#include "velocitybasis/sphericalbase.hpp"
#include "velocitybasis/sphericalharmonics.hpp"
#include <complex>

// ---- Linesource ----

MeltingCube::MeltingCube( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _sigmaS = settings->GetSigmaS(); }

MeltingCube::~MeltingCube() {}

VectorVector MeltingCube::GetScatteringXS( const Vector& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) );
}

VectorVector MeltingCube::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) ); }

std::vector<VectorVector> MeltingCube::GetExternalSource( const Vector& /*energies*/ ) {
    return std::vector<VectorVector>( 1u, std::vector<Vector>( _mesh->GetNumCells(), Vector( 1u, 0.0 ) ) );
}

// ---- MeltingCube_SN ----

MeltingCube_SN::MeltingCube_SN( Config* settings, Mesh* mesh ) : MeltingCube( settings, mesh ) {}

MeltingCube_SN::~MeltingCube_SN() {}

VectorVector MeltingCube_SN::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids          = _mesh->GetCellMidPoints();
    double kinetic_density = 0.0;
    double epsilon         = 1e-3;    // minimal value for first moment to avoid div by zero error
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];
        if( norm( cellMids[j] ) < 0.4 ) {
            kinetic_density = std::max( 0.2 * std::cos( norm( cellMids[j] ) * 4 ) * std::cos( norm( cellMids[j] ) * 4.0 ), epsilon );
        }
        else {
            kinetic_density = epsilon;
        }
        psi[j] = Vector( _settings->GetNQuadPoints(), kinetic_density );
    }
    return psi;
}

// ---- LineSource_PN ----

MeltingCube_Moment::MeltingCube_Moment( Config* settings, Mesh* mesh ) : MeltingCube( settings, mesh ) {}

MeltingCube_Moment::~MeltingCube_Moment() {}

VectorVector MeltingCube_Moment::SetupIC() {
    // Compute number of equations in the system

    // In case of PN, spherical basis is per default SPHERICAL_HARMONICS
    SphericalBase* tempBase  = SphericalBase::Create( _settings );
    unsigned ntotalEquations = tempBase->GetBasisSize();

    VectorVector psi( _mesh->GetNumCells(), Vector( ntotalEquations, 0 ) );    // zero could lead to problems?
    VectorVector cellMids = _mesh->GetCellMidPoints();

    Vector uIC( ntotalEquations, 0 );

    if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
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
    double epsilon = 1e-4;    // minimal value for first moment to avoid div by zero error
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double kinetic_density = 0.0;
        if( norm( cellMids[j] ) < 0.4 ) {
            kinetic_density = std::max( 0.2 * std::cos( norm( cellMids[j] ) * 4 ) * std::cos( norm( cellMids[j] ) * 4.0 ), epsilon );
        }
        else {
            kinetic_density = epsilon;
        }
        if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
            psi[j] = kinetic_density * uIC / uIC[0];    // Remember scaling
        }
        if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
            psi[j][0] = kinetic_density;
        }
    }
    delete tempBase;    // Only temporally needed
    return psi;
}

// ---- LineSource SN pseudo1D ----

MeltingCube_SN_1D::MeltingCube_SN_1D( Config* settings, Mesh* mesh ) : MeltingCube( settings, mesh ) {}

MeltingCube_SN_1D::~MeltingCube_SN_1D() {}

VectorVector MeltingCube_SN_1D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids          = _mesh->GetCellMidPoints();
    double epsilon         = 1e-3;    // minimal value for first moment to avoid div by zero error
    double kinetic_density = 0.0;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        if( x > -0.5 && x < 0 ) {
            kinetic_density = std::max( -x, epsilon );
        }
        else if( x > 0 && x < 0.5 ) {
            kinetic_density = std::max( -x, epsilon );
        }
        else {
            kinetic_density = epsilon;
        }
        psi[j] = kinetic_density;
    }
    return psi;
}

// ---- LineSource Moment pseudo1D ----

MeltingCube_Moment_1D::MeltingCube_Moment_1D( Config* settings, Mesh* mesh ) : MeltingCube( settings, mesh ) {}

MeltingCube_Moment_1D::~MeltingCube_Moment_1D() {}

VectorVector MeltingCube_Moment_1D::SetupIC() {
    double t       = 3.2e-4;    // pseudo time for gaussian smoothing (Approx to dirac impulse)
    double epsilon = 1e-3;      // minimal value for first moment to avoid div by zero error

    // In case of PN, spherical basis is per default SPHERICAL_HARMONICS
    if( _settings->GetSolverName() == PN_SOLVER || _settings->GetSolverName() == CSD_PN_SOLVER ) {
        // In case of PN, spherical basis is per default SPHERICAL_HARMONICS in 3 velocity dimensions
        SphericalHarmonics* tempBase = new SphericalHarmonics( _settings->GetMaxMomentDegree(), 3 );
        unsigned ntotalEquations     = tempBase->GetBasisSize();
        delete tempBase;
        VectorVector initialSolution( _mesh->GetNumCells(), Vector( ntotalEquations, 0.0 ) );    // zero could lead to problems?
        VectorVector cellMids  = _mesh->GetCellMidPoints();
        double kinetic_density = 0.0;
        for( unsigned idx_cell = 0; idx_cell < cellMids.size(); ++idx_cell ) {
            double x = cellMids[idx_cell][0];
            if( x > -0.4 && x < 0.4 ) {
                kinetic_density = std::max( 0.4 * std::cos( x * 4 ) * std::cos( x * 4.0 ), epsilon );
            }
            else {
                kinetic_density = epsilon;
            }
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

        if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
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
            double kinetic_density = 0.0;
            if( x > -0.4 && x < 0.4 ) {
                kinetic_density = std::max( 0.2 * std::cos( x * 4 ) * std::cos( x * 4.0 ), epsilon );
            }
            else {
                kinetic_density = epsilon;
            }
            if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
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
