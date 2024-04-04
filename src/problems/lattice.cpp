#include "problems/lattice.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "velocitybasis/sphericalbase.hpp"
#include "velocitybasis/sphericalharmonics.hpp"

// ---- Checkerboard Sn ----
// Constructor for Ckeckerboard case with Sn
Lattice_SN::Lattice_SN( Config* settings, Mesh* mesh, QuadratureBase* quad ) : ProblemBase( settings, mesh, quad ) {

    // Initialise scattering crosssections to 1 and absorption cross sections to 0
    _sigmaS = Vector( _mesh->GetNumCells(), 1. );
    _sigmaT = Vector( _mesh->GetNumCells(), 1. );

    // Initialize Quantities of interest
    _curAbsorptionLattice    = 0.0;
    _curMaxAbsorptionLattice = 0.0;
    _totalAbsorptionLattice  = 0.0;

    if( _settings->GetNLatticeAbsIndividual() == 49 && _settings->GetNLatticeScatterIndividual() == 49 ) {    // Individual values set
        auto log = spdlog::get( "event" );
        log->info( "| " );
        log->info( "| Lattice test case WITH individual scattering and absorption values for each block.  " );

        auto cellMids = _mesh->GetCellMidPoints();

        std::vector<double> scatteringValues = _settings->GetLatticeScatterIndividual();
        std::vector<double> absorptionValues = _settings->GetLatticeAbsorptionIndividual();

        for( unsigned j = 0; j < cellMids.size(); ++j ) {
            _sigmaS[j] = scatteringValues[GetBlockID( cellMids[j] )];
            _sigmaT[j] = absorptionValues[GetBlockID( cellMids[j] )] + scatteringValues[GetBlockID( cellMids[j] )];
        }
    }
    else {
        auto log = spdlog::get( "event" );
        log->info( "| " );
        log->info( "| Lattice test case WITHOUT individual scattering and absorption values for each block.  " );
        // For absorption cells: set scattering XS to 0 and absorption to 10
        auto cellMids = _mesh->GetCellMidPoints();
        for( unsigned j = 0; j < cellMids.size(); ++j ) {
            if( IsAbsorption( cellMids[j] ) ) {
                _sigmaS[j] = 0.0;
                _sigmaT[j] = _settings->GetLatticeAbsBlue();
            }
            else if( !IsSource( cellMids[j] ) ) {    // White block
                _sigmaS[j] = _settings->GetLatticeScatterWhite();
                _sigmaT[j] = _settings->GetLatticeScatterWhite();
            }
        }
    }

    SetGhostCells();
}

Lattice_SN::~Lattice_SN() {}

VectorVector Lattice_SN::GetScatteringXS( const Vector& /*energies */ ) { return VectorVector( 1u, _sigmaS ); }

VectorVector Lattice_SN::GetTotalXS( const Vector& /*energies */ ) { return VectorVector( 1u, _sigmaT ); }

std::vector<VectorVector> Lattice_SN::GetExternalSource( const Vector& /*energies*/ ) {
    VectorVector Q( _mesh->GetNumCells(), Vector( 1u, 0.0 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( IsSource( cellMids[j] ) ) Q[j] = _settings->GetSourceMagnitude();    // isotropic source
    }
    return std::vector<VectorVector>( 1u, Q );
}

VectorVector Lattice_SN::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-15 ) );
    return psi;
}

bool Lattice_SN::IsAbsorption( const Vector& pos ) const {
    // Check whether pos is inside absorbing squares
    double xy_corrector = -3.5;
    std::vector<double> lbounds{ 1 + xy_corrector, 2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector };
    std::vector<double> ubounds{ 2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector, 6 + xy_corrector };
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

bool Lattice_SN::IsAbsorption( double x, double y ) const {
    // Check whether pos is inside absorbing squares
    double xy_corrector = -3.5;
    std::vector<double> lbounds{ 1 + xy_corrector, 2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector };
    std::vector<double> ubounds{ 2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector, 6 + xy_corrector };
    for( unsigned k = 0; k < lbounds.size(); ++k ) {
        for( unsigned l = 0; l < lbounds.size(); ++l ) {
            if( ( l + k ) % 2 == 1 || ( k == 2 && l == 2 ) || ( k == 2 && l == 4 ) ) continue;
            if( x >= lbounds[k] && x <= ubounds[k] && y >= lbounds[l] && y <= ubounds[l] ) {
                return true;
            }
        }
    }
    return false;
}

unsigned Lattice_SN::GetBlockID( const Vector& pos ) const {
    double xy_corrector = 3.5;
    int block_x         = int( pos[0] + xy_corrector );
    int block_y         = int( pos[1] + xy_corrector );
    return (unsigned)( block_y * 7 + block_x );
}

bool Lattice_SN::IsSource( const Vector& pos ) const {
    // Check whether pos is part of source region
    if( pos[0] >= 3 - 3.5 && pos[0] <= 4 - 3.5 && pos[1] >= 3 - 3.5 && pos[1] <= 4 - 3.5 )
        return true;
    else
        return false;
}

void Lattice_SN::SetGhostCells() {
    // Loop over all cells. If its a Dirichlet boundary, add cell to dict with {cell_idx, boundary_value}
    auto cellBoundaries = _mesh->GetBoundaryTypes();
    std::map<int, Vector> ghostCellMap;

    QuadratureBase* quad = QuadratureBase::Create( _settings );
    VectorVector vq      = quad->GetPoints();
    unsigned nq          = quad->GetNq();

    Vector void_ghostcell( nq, 0.0 );

    for( unsigned idx_cell = 0; idx_cell < _mesh->GetNumCells(); idx_cell++ ) {
        if( cellBoundaries[idx_cell] == BOUNDARY_TYPE::NEUMANN || cellBoundaries[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) {
            ghostCellMap[idx_cell] = void_ghostcell;
        }
    }
    _ghostCells = ghostCellMap;

    delete quad;
}

const Vector& Lattice_SN::GetGhostCellValue( int idx_cell, const Vector& /* cell_sol */ ) { return _ghostCells[idx_cell]; }

// QOI getter
double Lattice_SN::GetCurAbsorptionLattice() { return _curAbsorptionLattice; }
double Lattice_SN::GetTotalAbsorptionLattice() { return _totalAbsorptionLattice; }
double Lattice_SN::GetMaxAbsorptionLattice() { return _curMaxAbsorptionLattice; }
// QOI setter
void Lattice_SN::ComputeTotalAbsorptionLattice( double dT ) { _totalAbsorptionLattice += _curAbsorptionLattice * dT; }

void Lattice_SN::ComputeCurrentAbsorptionLattice( const Vector& scalarFlux ) {
    _curAbsorptionLattice     = 0.0;
    unsigned nCells           = _mesh->GetNumCells();
    auto cellMids             = _mesh->GetCellMidPoints();
    std::vector<double> areas = _mesh->GetCellAreas();

    for( unsigned idx_cell = 0; idx_cell < nCells; idx_cell++ ) {
        if( IsAbsorption( cellMids[idx_cell] ) ) {
            _curAbsorptionLattice += scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * areas[idx_cell];
        }
    }
}
// TODO all absorption qois can be refactored in one function
void Lattice_SN::ComputeMaxAbsorptionLattice( const Vector& scalarFlux ) {
    unsigned nCells           = _mesh->GetNumCells();
    auto cellMids             = _mesh->GetCellMidPoints();
    std::vector<double> areas = _mesh->GetCellAreas();

    for( unsigned idx_cell = 0; idx_cell < nCells; idx_cell++ ) {
        if( IsAbsorption( cellMids[idx_cell] ) ) {
            if( _curMaxAbsorptionLattice < scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) )
                _curMaxAbsorptionLattice = scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] );
        }
    }
}

// ---- Checkerboard Moments ----

// Constructor for checkerboard case with Pn
Lattice_Moment::Lattice_Moment( Config* settings, Mesh* mesh, QuadratureBase* quad ) : ProblemBase( settings, mesh, quad ) {

    // Initialise crosssections = 1 (scattering)
    _sigmaS = Vector( _mesh->GetNumCells(), 1.0 );
    _sigmaT = Vector( _mesh->GetNumCells(), 1.0 );

    // for absorption regions change crosssections to all absorption
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isAbsorption( cellMids[j] ) ) {
            _sigmaS[j] = 0.0;
            _sigmaT[j] = 10.0;
        }
    }
}

Lattice_Moment::~Lattice_Moment() {}

VectorVector Lattice_Moment::GetScatteringXS( const Vector& /*energies*/ ) { return VectorVector( 1u, _sigmaS ); }

VectorVector Lattice_Moment::GetTotalXS( const Vector& /*energies*/ ) { return VectorVector( 1u, _sigmaT ); }

std::vector<VectorVector> Lattice_Moment::GetExternalSource( const Vector& /*energies*/ ) {
    // In case of PN, spherical basis is per default SPHERICAL_HARMONICS

    double integrationFactor = ( 4 * M_PI );
    if( _settings->GetDim() == 2 ) {
        integrationFactor = M_PI;
    }
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
                Q[j] = kinetic_density * uIC / uIC[0] / integrationFactor;    // Remember scaling
            }
            if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
                Q[j][0] = kinetic_density / integrationFactor;    // first bassis function is 1/ ( 4 * M_PI )
            }
        }
    }
    delete tempBase;    // Only temporally needed
    return std::vector<VectorVector>( 1u, Q );
}

VectorVector Lattice_Moment::SetupIC() {
    double integrationFactor = ( 4 * M_PI );
    if( _settings->GetDim() == 2 ) {
        integrationFactor = M_PI;
    }
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
    double kinetic_density = 1e-4;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
            initialSolution[j] = kinetic_density * tempIC / tempIC[0] / integrationFactor;    // Remember scaling
        }
        if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
            initialSolution[j][0] = kinetic_density / integrationFactor;    // first bassis function is 1/ ( 4 * M_PI )
        }
    }
    delete tempBase;    // Only temporally needed
    return initialSolution;
}

bool Lattice_Moment::isAbsorption( const Vector& pos ) const {
    // Check whether pos is in absorption region
    double xy_corrector = -3.5;
    std::vector<double> lbounds{ 1 + xy_corrector, 2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector };
    std::vector<double> ubounds{ 2 + xy_corrector, 3 + xy_corrector, 4 + xy_corrector, 5 + xy_corrector, 6 + xy_corrector };
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

bool Lattice_Moment::isSource( const Vector& pos ) const {
    // Check whether pos is in source region
    if( pos[0] >= 3 - 3.5 && pos[0] <= 4 - 3.5 && pos[1] >= 3 - 3.5 && pos[1] <= 4 - 3.5 )
        return true;
    else
        return false;
}
