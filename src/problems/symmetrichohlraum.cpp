#include "problems/symmetrichohlraum.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "quadratures/qgausslegendretensorized.hpp"
#include "quadratures/quadraturebase.hpp"
#include "solvers/csdpn_starmap_constants.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/interpolation.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "velocitybasis/sphericalbase.hpp"
#include "velocitybasis/sphericalharmonics.hpp"

SymmetricHohlraum::SymmetricHohlraum( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {
    _sigmaS = Vector( _mesh->GetNumCells(), 0.1 );    // white area default
    _sigmaT = Vector( _mesh->GetNumCells(), 0.1 );    // white area default

    _curAbsorptionHohlraumCenter       = 0.0;
    _curAbsorptionHohlraumVertical     = 0.0;
    _curAbsorptionHohlraumHorizontal   = 0.0;
    _totalAbsorptionHohlraumCenter     = 0.0;
    _totalAbsorptionHohlraumVertical   = 0.0;
    _totalAbsorptionHohlraumHorizontal = 0.0;
    _varAbsorptionHohlraumGreen        = 0.0;

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _mesh->GetNumCells(); idx_cell++ ) {
        // Assumption: Domain size is 1.3x1.3
        double x = _mesh->GetCellMidPoints()[idx_cell][0];
        double y = _mesh->GetCellMidPoints()[idx_cell][1];

        // red area left
        if( x < -0.6 && y > -0.4 && y < 0.4 ) {
            _sigmaS[idx_cell] = 95.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // red area right
        if( x > 0.6 && y > -0.4 && y < 0.4 ) {
            _sigmaS[idx_cell] = 95.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // green area 1 (lower boundary)
        if( x > -0.2 && x < -0.15 && y > -0.35 && y < 0.35 ) {
            _sigmaS[idx_cell] = 90.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // green area 2 (upper boundary)
        if( x > 0.15 && x < 0.2 && y > -0.35 && y < 0.35 ) {
            _sigmaS[idx_cell] = 90.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // green area 3 (left boundary)
        if( x > -0.2 && x < 0.2 && y > -0.4 && y < -0.35 ) {
            _sigmaS[idx_cell] = 90.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // green area 4 (right boundary)
        if( x > -0.2 && x < 0.2 && y > 0.35 && y < 0.4 ) {
            _sigmaS[idx_cell] = 90.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // blue checkered area
        if( x > -0.15 && x < 0.15 && y > -0.35 && y < 0.35 ) {
            _sigmaS[idx_cell] = 50.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // black area (upper and lower boundary)
        if( y > 0.6 || y < -0.6 ) {
            _sigmaS[idx_cell] = 0.0;
            _sigmaT[idx_cell] = 100.0;
        }
    }
    SetGhostCells();
}

SymmetricHohlraum::~SymmetricHohlraum() {}

std::vector<VectorVector> SymmetricHohlraum::GetExternalSource( const Vector& energies ) {
    VectorVector Q( _mesh->GetNumCells(), Vector( 1u, 0.0 ) );
    return std::vector<VectorVector>( 1u, Q );
}

VectorVector SymmetricHohlraum::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    VectorVector cellMids = _mesh->GetCellMidPoints();

    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        psi[j] = 1e-4;    // zero initial condition
    }
    return psi;
}

void SymmetricHohlraum::SetGhostCells() {
    // Loop over all cells. If its a Dirichlet boundary, add cell to dict with {cell_idx, boundary_value}
    auto cellBoundaries = _mesh->GetBoundaryTypes();
    std::map<int, Vector> ghostCellMap;

    QuadratureBase* quad = QuadratureBase::Create( _settings );
    VectorVector vq      = quad->GetPoints();
    unsigned nq          = quad->GetNq();

    Vector left_inflow( nq, 0.0 );
    Vector right_inflow( nq, 0.0 );
    Vector vertical_flow( nq, 0.0 );

    for( unsigned idx_q = 0; idx_q < nq; idx_q++ ) {
        if( vq[idx_q][0] > 0.0 ) left_inflow[idx_q] = 1.0;
        if( vq[idx_q][0] < 0.0 ) right_inflow[idx_q] = 1.0;
    }

    for( unsigned idx_cell = 0; idx_cell < _mesh->GetNumCells(); idx_cell++ ) {
        double x = _mesh->GetCellMidPoints()[idx_cell][0];
        double y = _mesh->GetCellMidPoints()[idx_cell][1];

        if( cellBoundaries[idx_cell] == BOUNDARY_TYPE::NEUMANN || cellBoundaries[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) {
            if( y < -0.6 )
                ghostCellMap.insert( { idx_cell, vertical_flow } );
            else if( y > 0.6 )
                ghostCellMap.insert( { idx_cell, vertical_flow } );
            else if( x < -0.6 )
                ghostCellMap.insert( { idx_cell, left_inflow } );
            else if( x > 0.6 )
                ghostCellMap.insert( { idx_cell, right_inflow } );
        }
    }
    _ghostCells = ghostCellMap;

    delete quad;
}

const Vector& SymmetricHohlraum::GetGhostCellValue( int idx_cell, const Vector& cell_sol ) { return _ghostCells[idx_cell]; }

VectorVector SymmetricHohlraum::GetScatteringXS( const Vector& energies ) { return VectorVector( 1u, _sigmaS ); }

VectorVector SymmetricHohlraum::GetTotalXS( const Vector& energies ) { return VectorVector( 1u, _sigmaT ); }

void SymmetricHohlraum::ComputeCurrentAbsorptionHohlraum( const Vector& scalarFlux ) {
    _curAbsorptionHohlraumCenter     = 0.0;    // Green and blue areas of symmetric hohlraum
    _curAbsorptionHohlraumVertical   = 0.0;    // Red areas of symmetric hohlraum
    _curAbsorptionHohlraumHorizontal = 0.0;    // Black areas of symmetric hohlraum

    unsigned nCells           = _mesh->GetNumCells();
    auto cellMids             = _mesh->GetCellMidPoints();
    std::vector<double> areas = _mesh->GetCellAreas();

    for( unsigned idx_cell = 0; idx_cell < nCells; idx_cell++ ) {
        double x = _mesh->GetCellMidPoints()[idx_cell][0];
        double y = _mesh->GetCellMidPoints()[idx_cell][1];

        if( x > -0.2 && x < 0.2 && y > -0.35 && y < 0.35 ) {
            _curAbsorptionHohlraumCenter += scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * areas[idx_cell];
        }
        if( ( x < -0.6 && y > -0.4 && y < 0.4 ) || ( x > 0.6 && y > -0.4 && y < 0.4 ) ) {
            _curAbsorptionHohlraumVertical += scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * areas[idx_cell];
        }
        if( y > 0.6 || y < -0.6 ) {
            _curAbsorptionHohlraumHorizontal += scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * areas[idx_cell];
        }
    }
}

void SymmetricHohlraum::ComputeTotalAbsorptionHohlraum( double dT ) {
    _totalAbsorptionHohlraumCenter += _curAbsorptionHohlraumCenter * dT;
    _totalAbsorptionHohlraumVertical += _curAbsorptionHohlraumVertical * dT;
    _totalAbsorptionHohlraumHorizontal += _curAbsorptionHohlraumHorizontal * dT;
}

void SymmetricHohlraum::ComputeVarAbsorptionGreen( const Vector& scalarFlux ) {
    bool green1, green2, green3, green4;
    double a_g                  = 0.0;
    _varAbsorptionHohlraumGreen = 0.0;
    double x, y;

    unsigned nCells           = _mesh->GetNumCells();
    auto cellMids             = _mesh->GetCellMidPoints();
    std::vector<double> areas = _mesh->GetCellAreas();

    for( unsigned idx_cell = 0; idx_cell < nCells; ++idx_cell ) {
        x = _mesh->GetCellMidPoints()[idx_cell][0];
        y = _mesh->GetCellMidPoints()[idx_cell][1];

        green1 = x > -0.2 && x < -0.15 && y > -0.35 && y < 0.35;    // green area 1 (lower boundary)
        green2 = x > 0.15 && x < 0.2 && y > -0.35 && y < 0.35;      // green area 2 (upper boundary)
        green3 = x > -0.2 && x < 0.2 && y > -0.4 && y < -0.35;      // green area 3 (left boundary)
        green4 = x > -0.2 && x < 0.2 && y > 0.35 && y < 0.4;        // green area 4 (right boundary)

        if( green1 || green2 || green3 || green4 ) {
            a_g += ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * scalarFlux[idx_cell] * areas[idx_cell];
        }
    }
    for( unsigned idx_cell = 0; idx_cell < nCells; ++idx_cell ) {
        x = _mesh->GetCellMidPoints()[idx_cell][0];
        y = _mesh->GetCellMidPoints()[idx_cell][1];

        green1 = x > -0.2 && x < -0.15 && y > -0.35 && y < 0.35;    // green area 1 (lower boundary)
        green2 = x > 0.15 && x < 0.2 && y > -0.35 && y < 0.35;      // green area 2 (upper boundary)
        green3 = x > -0.2 && x < 0.2 && y > -0.4 && y < -0.35;      // green area 3 (left boundary)
        green4 = x > -0.2 && x < 0.2 && y > 0.35 && y < 0.4;        // green area 4 (right boundary)

        if( green1 || green2 || green3 || green4 ) {
            _varAbsorptionHohlraumGreen += ( a_g - scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) ) *
                                           ( a_g - scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) ) * areas[idx_cell];
        }
    }
}

// -------------- Moment Symmetric Hohlraum ---------------

SymmetricHohlraum_Moment::SymmetricHohlraum_Moment( Config* settings, Mesh* mesh ) : SymmetricHohlraum( settings, mesh ) {}

SymmetricHohlraum_Moment::~SymmetricHohlraum_Moment() {}

std::vector<VectorVector> SymmetricHohlraum_Moment::GetExternalSource( const Vector& energies ) {
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
    double kinetic_density = 0.0;    //_settings->GetSourceMagnitude();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( cellMids[j][0] < 0.05 ) {
            if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS || _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS_ROTATED ) {
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

VectorVector SymmetricHohlraum_Moment::SetupIC() {
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
            tempIC += w[idx_quad] * moments[idx_quad];
        }
        delete quad;
    }
    // Initial condition is dirac impulse at (x,y) = (0,0) ==> constant in angle ==> all moments - exept first - are zero.
    double kinetic_density = 1e-4;
    // std::vector<BOUNDARY_TYPE> _boundaryCells;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {

        // boundary condition: Source on left side
        if( cellMids[j][0] < 0.0 && ( cellMids[j][1] > 0.0 && cellMids[j][1] < 1.3 ) ) {    // test case uses ghost cells
            kinetic_density = _settings->GetSourceMagnitude();
            _mesh->SetBoundaryType( j, DIRICHLET );
        }
        else {
            kinetic_density = 1e-4;
        }

        if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS || _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS_ROTATED ) {
            initialSolution[j] = kinetic_density * tempIC / tempIC[0] / integrationFactor;    // Remember scaling
        }
        if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
            initialSolution[j][0] = kinetic_density / integrationFactor;    // first bassis function is 1/ ( 4 * M_PI )
        }
    }
    delete tempBase;    // Only temporally needed
    return initialSolution;
}
