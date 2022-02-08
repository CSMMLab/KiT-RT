#include "problems/checkerboard.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/sphericalbase.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

// ---- Checkerboard Sn ----
// Constructor for Ckeckerboard case with Sn
Checkerboard_SN::Checkerboard_SN( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {

    // Initialise crosssections to 1
    _scatteringXS = Vector( _mesh->GetNumCells(), 1.0 );
    _totalXS      = Vector( _mesh->GetNumCells(), 1.0 );

    // For absorption cells: set scattering XS to 0 and absorption to 10
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isAbsorption( cellMids[j] ) ) {
            _scatteringXS[j] = 0.0;
            _totalXS[j]      = 10.0;
        }
    }
}

Checkerboard_SN::~Checkerboard_SN() {}

VectorVector Checkerboard_SN::GetScatteringXS( const Vector& energies ) { return VectorVector( energies.size(), _scatteringXS ); }

VectorVector Checkerboard_SN::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), _totalXS ); }

std::vector<VectorVector> Checkerboard_SN::GetExternalSource( const Vector& /*energies*/ ) {
    VectorVector Q( _mesh->GetNumCells(), Vector( 1u, 0.0 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isSource( cellMids[j] ) ) Q[j] = _settings->GetSourceMagnitude() / ( 4 * M_PI );    // isotropic source
    }
    return std::vector<VectorVector>( 1u, Q );
}

VectorVector Checkerboard_SN::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    return psi;
}

bool Checkerboard_SN::isAbsorption( const Vector& pos ) const {
    // Check whether pos is inside absorbing squares
    std::vector<double> lbounds{ 1, 2, 3, 4, 5 };
    std::vector<double> ubounds{ 2, 3, 4, 5, 6 };
    for( unsigned k = 0; k < lbounds.size(); ++k ) {
        for( unsigned l = 0; l < lbounds.size(); ++l ) {
            if( ( l + k ) % 2 == 1 || ( k == 2 && l == 2 ) || ( k == 2 && l == 4 ) ) continue;
            if( pos[0] >= lbounds[k] && pos[0] <= ubounds[k] && pos[1] >= lbounds[l] && pos[1] <= ubounds[l] ) {
                return true;
            }
        }
    }
    return false;
}

bool Checkerboard_SN::isSource( const Vector& pos ) const {
    // Check whether pos is part of source region
    if( pos[0] >= 3 && pos[0] <= 4 && pos[1] >= 3 && pos[1] <= 4 )
        return true;
    else
        return false;
}

// ---- Checkerboard Moments ----

// Constructor for checkerboard case with Pn
Checkerboard_Moment::Checkerboard_Moment( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {

    // Initialise crosssections = 1 (scattering)
    _scatteringXS = Vector( _mesh->GetNumCells(), 1.0 );
    _totalXS      = Vector( _mesh->GetNumCells(), 1.0 );

    // for absorption regions change crosssections to all absorption
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isAbsorption( cellMids[j] ) ) {
            _scatteringXS[j] = 0.0;
            _totalXS[j]      = 10.0;
        }
    }
}

Checkerboard_Moment::~Checkerboard_Moment() {}

VectorVector Checkerboard_Moment::GetScatteringXS( const Vector& energies ) { return VectorVector( energies.size(), _scatteringXS ); }

VectorVector Checkerboard_Moment::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), _totalXS ); }

std::vector<VectorVector> Checkerboard_Moment::GetExternalSource( const Vector& /*energies*/ ) {
    // In case of PN, spherical basis is per default SPHERICAL_HARMONICS
    SphericalBase* tempBase  = SphericalBase::Create( _settings );
    unsigned ntotalEquations = tempBase->GetBasisSize();

    VectorVector Q( _mesh->GetNumCells(), Vector( ntotalEquations, 0.0 ) );    // zero could lead to problems?
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
    double kinetic_density = _settings->GetSourceMagnitude();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isSource( cellMids[j] ) ) {
            if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
                Q[j] = kinetic_density * uIC / uIC[0];    // Remember scaling
            }
            if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
                Q[j][0] = kinetic_density / std::sqrt( 4 * M_PI );
            }
        }
    }
    delete tempBase;    // Only temporally needed
    return std::vector<VectorVector>( 1u, Q );
}

VectorVector Checkerboard_Moment::SetupIC() {
    // In case of PN, spherical basis is per default SPHERICAL_HARMONICS
    SphericalBase* tempBase  = SphericalBase::Create( _settings );
    unsigned ntotalEquations = tempBase->GetBasisSize();

    VectorVector initialSolution( _mesh->GetNumCells(), Vector( ntotalEquations, 0 ) );    // zero could lead to problems?
    VectorVector cellMids = _mesh->GetCellMidPoints();

    Vector tempIC( ntotalEquations, 0 );

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
            tempIC += w[idx_quad] * moments[idx_quad];
        }
        delete quad;
    }
    // Initial condition is dirac impulse at (x,y) = (0,0) ==> constant in angle ==> all moments - exept first - are zero.
    double kinetic_density = 1e-3;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
            initialSolution[j] = kinetic_density * tempIC / tempIC[0];    // Remember scaling
        }
        if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
            initialSolution[j][0] = kinetic_density;
        }
    }
    delete tempBase;    // Only temporally needed
    return initialSolution;
}

bool Checkerboard_Moment::isAbsorption( const Vector& pos ) const {
    // Check whether pos is in absorption region
    std::vector<double> lbounds{ 1, 2, 3, 4, 5 };
    std::vector<double> ubounds{ 2, 3, 4, 5, 6 };
    for( unsigned k = 0; k < lbounds.size(); ++k ) {
        for( unsigned l = 0; l < lbounds.size(); ++l ) {
            if( ( l + k ) % 2 == 1 || ( k == 2 && l == 2 ) || ( k == 2 && l == 4 ) ) continue;
            if( pos[0] >= lbounds[k] && pos[0] <= ubounds[k] && pos[1] >= lbounds[l] && pos[1] <= ubounds[l] ) {
                return true;
            }
        }
    }
    return false;
}

bool Checkerboard_Moment::isSource( const Vector& pos ) const {
    // Check whether pos is in source region
    if( pos[0] >= 3 && pos[0] <= 4 && pos[1] >= 3 && pos[1] <= 4 )
        return true;
    else
        return false;
}

// ---- Checkerboard SN 1D ----
// Constructor for Ckeckerboard case with Sn
Checkerboard_SN_1D::Checkerboard_SN_1D( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {

    // Initialise crosssections to 1
    _scatteringXS = Vector( _mesh->GetNumCells(), 1.0 );
    _totalXS      = Vector( _mesh->GetNumCells(), 1.0 );

    // For absorption cells: set scattering XS to 0 and absorption to 10
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isAbsorption( cellMids[j] ) ) {
            _scatteringXS[j] = 0.0;
            _totalXS[j]      = 10.0;
        }
    }
}

Checkerboard_SN_1D::~Checkerboard_SN_1D() {}

VectorVector Checkerboard_SN_1D::GetScatteringXS( const Vector& energies ) { return VectorVector( energies.size(), _scatteringXS ); }

VectorVector Checkerboard_SN_1D::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), _totalXS ); }

std::vector<VectorVector> Checkerboard_SN_1D::GetExternalSource( const Vector& /*energies*/ ) {
    VectorVector Q( _mesh->GetNumCells(), Vector( 1u, 0.0 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isSource( cellMids[j] ) ) Q[j] = _settings->GetSourceMagnitude() / ( 4 * M_PI );    // isotropic source
    }
    return std::vector<VectorVector>( 1u, Q );
}

VectorVector Checkerboard_SN_1D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    return psi;
}

bool Checkerboard_SN_1D::isAbsorption( const Vector& pos ) const {
    // Check whether pos is inside absorbing squares
    // domain from 0 to 7, absorption block is between 1 and 2
    if( ( pos[0] >= 1 && pos[0] <= 2 ) || ( pos[0] >= 6.5 && pos[0] <= 7 ) ) {
        return true;
    }
    return false;
}

bool Checkerboard_SN_1D::isSource( const Vector& pos ) const {
    // Check whether pos is part of source region
    if( pos[0] >= 3 && pos[0] <= 4 )
        return true;
    else
        return false;
}

// --- Moment version 1d ---

Checkerboard_Moment_1D::Checkerboard_Moment_1D( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {

    // Initialise crosssections to 1
    _scatteringXS = Vector( _mesh->GetNumCells(), 1.0 );
    _totalXS      = Vector( _mesh->GetNumCells(), 1.0 );

    // For absorption cells: set scattering XS to 0 and absorption to 10
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isAbsorption( cellMids[j] ) ) {
            _scatteringXS[j] = 0.0;
            _totalXS[j]      = 10.0;
        }
    }
}

Checkerboard_Moment_1D::~Checkerboard_Moment_1D() {}

VectorVector Checkerboard_Moment_1D::GetScatteringXS( const Vector& energies ) { return VectorVector( energies.size(), _scatteringXS ); }

VectorVector Checkerboard_Moment_1D::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), _totalXS ); }

std::vector<VectorVector> Checkerboard_Moment_1D::GetExternalSource( const Vector& /*energies*/ ) {
    if( _settings->GetSolverName() == PN_SOLVER || _settings->GetSolverName() == CSD_PN_SOLVER ) {
        // In case of PN, spherical basis is per default SPHERICAL_HARMONICS in 3 velocity dimensions
        SphericalBase* tempBase  = new SphericalHarmonics( _settings->GetMaxMomentDegree(), 3 );
        unsigned ntotalEquations = tempBase->GetBasisSize();
        delete tempBase;
        VectorVector Q( _mesh->GetNumCells(), Vector( ntotalEquations, 0.0 ) );
        double kinetic_density = _settings->GetSourceMagnitude();
        VectorVector cellMids  = _mesh->GetCellMidPoints();
        for( unsigned idx_cell = 0; idx_cell < cellMids.size(); ++idx_cell ) {
            if( isSource( cellMids[idx_cell] ) ) {
                Q[idx_cell][0] = kinetic_density / std::sqrt( 4 * M_PI );
            }
        }
        return std::vector<VectorVector>( 1u, Q );
    }
    else {
        // In case of PN, spherical basis is per default SPHERICAL_HARMONICS
        SphericalBase* tempBase  = SphericalBase::Create( _settings );
        unsigned ntotalEquations = tempBase->GetBasisSize();

        VectorVector Q( _mesh->GetNumCells(), Vector( ntotalEquations, 0.0 ) );    // zero could lead to problems?
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
        double kinetic_density = _settings->GetSourceMagnitude();
        for( unsigned j = 0; j < cellMids.size(); ++j ) {
            if( isSource( cellMids[j] ) ) {
                if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
                    Q[j] = kinetic_density * uIC / uIC[0];    // Remember scaling
                }
                if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
                    Q[j][0] = kinetic_density / std::sqrt( 4 * M_PI );
                }
            }
        }
        delete tempBase;    // Only temporally needed
        return std::vector<VectorVector>( 1u, Q );
    }
}

VectorVector Checkerboard_Moment_1D::SetupIC() {
    if( _settings->GetSolverName() == PN_SOLVER || _settings->GetSolverName() == CSD_PN_SOLVER ) {
        // In case of PN, spherical basis is per default SPHERICAL_HARMONICS in 3 velocity dimensions
        SphericalBase* tempBase  = new SphericalHarmonics( _settings->GetMaxMomentDegree(), 3 );
        unsigned ntotalEquations = tempBase->GetBasisSize();
        delete tempBase;
        double epsilon = 1e-3;
        VectorVector initialSolution( _mesh->GetNumCells(), Vector( ntotalEquations, 0.0 ) );    // zero could lead to problems?
        VectorVector cellMids = _mesh->GetCellMidPoints();
        for( unsigned idx_cell = 0; idx_cell < cellMids.size(); ++idx_cell ) {
            initialSolution[idx_cell][0] = epsilon;
        }
        return initialSolution;
    }
    else {
        SphericalBase* tempBase  = SphericalBase::Create( _settings );
        unsigned ntotalEquations = tempBase->GetBasisSize();

        VectorVector initialSolution( _mesh->GetNumCells(), Vector( ntotalEquations, 0 ) );    // zero could lead to problems?
        VectorVector cellMids = _mesh->GetCellMidPoints();
        Vector tempIC( ntotalEquations, 0 );
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
                tempIC += w[idx_quad] * moments[idx_quad];
            }
            delete quad;
        }
        // Initial condition is dirac impulse at (x,y) = (0,0) ==> constant in angle ==> all moments - exept first - are zero.
        double kinetic_density = 1e-4;
        for( unsigned j = 0; j < cellMids.size(); ++j ) {
            if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
                initialSolution[j] = kinetic_density * tempIC / tempIC[0];    // Remember scaling
            }
            if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
                initialSolution[j][0] = kinetic_density;
            }
        }
        delete tempBase;    // Only temporally needed
        return initialSolution;
    }
}

bool Checkerboard_Moment_1D::isAbsorption( const Vector& pos ) const {
    // Check whether pos is inside absorbing squares
    // domain from 0 to 7, absorption block is between 1 and 2
    if( ( pos[0] >= 1 && pos[0] <= 2 ) || ( pos[0] >= 6.5 && pos[0] <= 7 ) ) {
        return true;
    }
    return false;
}

bool Checkerboard_Moment_1D::isSource( const Vector& pos ) const {
    // Check whether pos is part of source region
    if( pos[0] >= 3 && pos[0] <= 4 )
        return true;
    else
        return false;
}
