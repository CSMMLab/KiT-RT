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

SymmetricHohlraum::SymmetricHohlraum( Config* settings, Mesh* mesh, QuadratureBase* quad ) : ProblemBase( settings, mesh, quad ) {
    _sigmaS = Vector( _mesh->GetNumCells(), 0.1 );    // white area default
    _sigmaT = Vector( _mesh->GetNumCells(), 0.1 );    // white area default

    // Design parameters defining the hohlraum geometry

    // Red
    _redLeftTop        = _settings->GetPosRedLeftTopHohlraum();
    _redLeftBottom     = _settings->GetPosRedLeftBottomHohlraum();
    _redRightTop       = _settings->GetPosRedRightTopHohlraum();
    _redRightBottom    = _settings->GetPosRedRightBottomHohlraum();
    _posRedLeftBorder  = _settings->GetPosRedLeftBorderHohlraum();
    _posRedRightBorder = _settings->GetPosRedRightBorderHohlraum();

    // Green
    _widthGreen     = 0.4;
    _heightGreen    = 0.8;
    _thicknessGreen = 0.05;
    _centerGreen    = { _settings->GetPosXCenterGreenHohlraum(), _settings->GetPosYCenterGreenHohlraum() };

    _cornerUpperLeftGreen  = { _centerGreen[0] - _widthGreen / 2.0, _centerGreen[1] + _heightGreen / 2.0 };
    _cornerLowerLeftGreen  = { _centerGreen[0] - _widthGreen / 2.0, _centerGreen[1] - _heightGreen / 2.0 };
    _cornerUpperRightGreen = { _centerGreen[0] + _widthGreen / 2.0, _centerGreen[1] + _heightGreen / 2.0 };
    _cornerLowerRightGreen = { _centerGreen[0] + _widthGreen / 2.0, _centerGreen[1] - _heightGreen / 2.0 };

    // QOIS
    _curAbsorptionHohlraumCenter       = 0.0;
    _curAbsorptionHohlraumVertical     = 0.0;
    _curAbsorptionHohlraumHorizontal   = 0.0;
    _totalAbsorptionHohlraumCenter     = 0.0;
    _totalAbsorptionHohlraumVertical   = 0.0;
    _totalAbsorptionHohlraumHorizontal = 0.0;
    _varAbsorptionHohlraumGreen        = 0.0;

    _probingCells = {
        _mesh->GetCellOfKoordinate( -0.4, 0. ),
        _mesh->GetCellOfKoordinate( 0.4, 0. ),
        _mesh->GetCellOfKoordinate( 0., -0.5 ),
        _mesh->GetCellOfKoordinate( 0., 0.5 ),
    };

    _probingMoments         = VectorVector( 4, Vector( 3, 0.0 ) );
    _nProbingCellsLineGreen = _settings->GetNumProbingCellsLineHohlraum();

    SetProbingCellsLineGreen();
    _absorptionValsIntegrated    = std::vector<double>( _nProbingCellsLineGreen, 0.0 );
    _varAbsorptionValsIntegrated = std::vector<double>( _nProbingCellsLineGreen, 0.0 );

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _mesh->GetNumCells(); idx_cell++ ) {
        // Assumption: Domain size is 1.3x1.3
        double x = _mesh->GetCellMidPoints()[idx_cell][0];
        double y = _mesh->GetCellMidPoints()[idx_cell][1];

        // red area left
        if( x < _posRedLeftBorder && y > _redLeftBottom && y < _redLeftTop ) {
            _sigmaS[idx_cell] = 95.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // red area right
        if( x > _posRedRightBorder && y > _redRightBottom && y < _redRightTop ) {
            _sigmaS[idx_cell] = 95.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // green and blue area
        if( x > -0.2 + _centerGreen[0] && x < 0.2 + _centerGreen[0] && y > -0.4 + _centerGreen[1] && y < 0.4 + _centerGreen[1] ) {
            _sigmaS[idx_cell] = 90.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // blue checkered area (overwrites part of green n blue area)
        if( x > -0.15 + _centerGreen[0] && x < 0.15 + _centerGreen[0] && y > -0.35 + _centerGreen[1] && y < 0.35 + _centerGreen[1] ) {
            _sigmaS[idx_cell] = 50.0;
            _sigmaT[idx_cell] = 100.0;
        }
        // black area (upper and lower boundary)
        if( y > 0.6 || y < -0.6 ) {
            _sigmaS[idx_cell] = 100.0;
            _sigmaT[idx_cell] = 100.0;
        }
    }
    SetGhostCells();
}

SymmetricHohlraum::~SymmetricHohlraum() {}

std::vector<VectorVector> SymmetricHohlraum::GetExternalSource( const Vector& /* energies */ ) {
    VectorVector Q( _mesh->GetNumCells(), Vector( 1u, 0.0 ) );
    return std::vector<VectorVector>( 1u, Q );
}

VectorVector SymmetricHohlraum::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) );
    VectorVector cellMids = _mesh->GetCellMidPoints();

    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        psi[j] = 0.0;    // zero initial condition
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

const Vector& SymmetricHohlraum::GetGhostCellValue( int idx_cell, const Vector& /* cell_sol */ ) { return _ghostCells[idx_cell]; }

VectorVector SymmetricHohlraum::GetScatteringXS( const Vector& /* energies */ ) { return VectorVector( 1u, _sigmaS ); }

VectorVector SymmetricHohlraum::GetTotalXS( const Vector& /* energies */ ) { return VectorVector( 1u, _sigmaT ); }

void SymmetricHohlraum::ComputeCurrentAbsorptionHohlraum( const Vector& scalarFlux ) {
    _curAbsorptionHohlraumCenter     = 0.0;    // Green and blue areas of symmetric hohlraum
    _curAbsorptionHohlraumVertical   = 0.0;    // Red areas of symmetric hohlraum
    _curAbsorptionHohlraumHorizontal = 0.0;    // Black areas of symmetric hohlraum

    unsigned nCells           = _mesh->GetNumCells();
    auto cellMids             = _mesh->GetCellMidPoints();
    std::vector<double> areas = _mesh->GetCellAreas();

#pragma omp parallel for default( shared )                                                                                                           \
    reduction( + : _curAbsorptionHohlraumCenter, _curAbsorptionHohlraumVertical, _curAbsorptionHohlraumHorizontal )
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
    double a_g                  = 0.0;
    _varAbsorptionHohlraumGreen = 0.0;

    unsigned nCells           = _mesh->GetNumCells();
    auto cellMids             = _mesh->GetCellMidPoints();
    std::vector<double> areas = _mesh->GetCellAreas();

#pragma omp parallel for default( shared ) reduction( + : a_g )
    for( unsigned idx_cell = 0; idx_cell < nCells; ++idx_cell ) {
        double x = cellMids[idx_cell][0];
        double y = cellMids[idx_cell][1];
        // green area 1 (lower boundary)
        bool green1 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < -0.15 + _settings->GetPosXCenterGreenHohlraum() &&
                      y > -0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.35 + _settings->GetPosYCenterGreenHohlraum();
        // green area 2 (upper boundary)
        bool green2 = x > 0.15 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                      y > -0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.35 + _settings->GetPosYCenterGreenHohlraum();
        // green area 3 (left boundary)
        bool green3 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                      y > -0.4 + _settings->GetPosYCenterGreenHohlraum() && y < -0.35 + _settings->GetPosYCenterGreenHohlraum();
        // green area 4 (right boundary)
        bool green4 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                      y > 0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.4 + _settings->GetPosYCenterGreenHohlraum();

        if( green1 || green2 || green3 || green4 ) {
            a_g += ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) * scalarFlux[idx_cell] * areas[idx_cell];
        }
    }

#pragma omp parallel for default( shared ) reduction( + : _varAbsorptionHohlraumGreen )
    for( unsigned idx_cell = 0; idx_cell < nCells; ++idx_cell ) {
        double x = cellMids[idx_cell][0];
        double y = cellMids[idx_cell][1];
        // green area 1 (lower boundary)
        bool green1 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < -0.15 + _settings->GetPosXCenterGreenHohlraum() &&
                      y > -0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.35 + _settings->GetPosYCenterGreenHohlraum();
        // green area 2 (upper boundary)
        bool green2 = x > 0.15 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                      y > -0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.35 + _settings->GetPosYCenterGreenHohlraum();
        // green area 3 (left boundary)
        bool green3 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                      y > -0.4 + _settings->GetPosYCenterGreenHohlraum() && y < -0.35 + _settings->GetPosYCenterGreenHohlraum();
        // green area 4 (right boundary)
        bool green4 = x > -0.2 + _settings->GetPosXCenterGreenHohlraum() && x < 0.2 + _settings->GetPosXCenterGreenHohlraum() &&
                      y > 0.35 + _settings->GetPosYCenterGreenHohlraum() && y < 0.4 + _settings->GetPosYCenterGreenHohlraum();

        if( green1 || green2 || green3 || green4 ) {
            _varAbsorptionHohlraumGreen += ( a_g - scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) ) *
                                           ( a_g - scalarFlux[idx_cell] * ( _sigmaT[idx_cell] - _sigmaS[idx_cell] ) ) * areas[idx_cell];
        }
    }
}

void SymmetricHohlraum::ComputeCurrentProbeMoment( const VectorVector& solution ) {
    const VectorVector& quadPoints = _quad->GetPoints();
    const Vector& weights          = _quad->GetWeights();
    unsigned nq                    = _quad->GetNq();

    for( unsigned idx_cell = 0; idx_cell < 4; idx_cell++ ) {    // Loop over probing cells
        _probingMoments[idx_cell][0] = blaze::dot( solution[_probingCells[idx_cell]], weights );
        _probingMoments[idx_cell][1] = 0.0;
        _probingMoments[idx_cell][2] = 0.0;

        for( unsigned idx_quad = 0; idx_quad < nq; idx_quad++ ) {
            _probingMoments[idx_cell][1] += quadPoints[idx_quad][0] * solution[_probingCells[idx_cell]][idx_quad] * weights[idx_quad];
            _probingMoments[idx_cell][2] += quadPoints[idx_quad][2] * solution[_probingCells[idx_cell]][idx_quad] * weights[idx_quad];
        }
    }
}

void SymmetricHohlraum::SetProbingCellsLineGreen() {

    double verticalLineWidth   = std::abs( _cornerUpperLeftGreen[1] - _cornerLowerLeftGreen[1] - _thicknessGreen );
    double horizontalLineWidth = std::abs( _cornerUpperLeftGreen[0] - _cornerUpperRightGreen[0] );

    // double dx = 2 * ( horizontalLineWidth + verticalLineWidth ) / ( (double)_nProbingCellsLineGreen );

    unsigned nHorizontalProbingCells =
        (unsigned)std::ceil( _nProbingCellsLineGreen / 2 * ( horizontalLineWidth / ( horizontalLineWidth + verticalLineWidth ) ) );
    unsigned nVerticalProbingCells = _nProbingCellsLineGreen - nHorizontalProbingCells;

    _probingCellsLineGreen = std::vector<unsigned>( _nProbingCellsLineGreen );

    std::vector<double> p1 = { _cornerUpperLeftGreen[0] + _thicknessGreen / 2.0, _cornerUpperLeftGreen[1] - _thicknessGreen / 2.0 };
    std::vector<double> p2 = { _cornerLowerLeftGreen[0] + _thicknessGreen / 2.0, _cornerLowerLeftGreen[1] + _thicknessGreen / 2.0 };
    std::vector<double> p3 = { _cornerUpperRightGreen[0] - _thicknessGreen / 2.0, _cornerUpperRightGreen[1] - _thicknessGreen / 2.0 };
    std::vector<double> p4 = { _cornerLowerRightGreen[0] - _thicknessGreen / 2.0, _cornerLowerRightGreen[1] + _thicknessGreen / 2.0 };

    // Sample points on each side of the rectangle
    std::vector<unsigned> side1 = linspace2D( p1, p2, nVerticalProbingCells );
    std::vector<unsigned> side2 = linspace2D( p2, p3, nHorizontalProbingCells );
    std::vector<unsigned> side3 = linspace2D( p3, p4, nVerticalProbingCells );
    std::vector<unsigned> side4 = linspace2D( p4, p1, nHorizontalProbingCells );

    // printf( "here" );
    //  Combine the points from each side
    _probingCellsLineGreen.insert( _probingCellsLineGreen.end(), side1.begin(), side1.end() );
    _probingCellsLineGreen.insert( _probingCellsLineGreen.end(), side2.begin(), side2.end() );
    _probingCellsLineGreen.insert( _probingCellsLineGreen.end(), side3.begin(), side3.end() );
    _probingCellsLineGreen.insert( _probingCellsLineGreen.end(), side4.begin(), side4.end() );
}

void SymmetricHohlraum::ComputeQOIsGreenProbingLine( const Vector& scalarFlux ) {

    double verticalLineWidth   = std::abs( _cornerUpperLeftGreen[1] - _cornerLowerLeftGreen[1] - _thicknessGreen );
    double horizontalLineWidth = std::abs( _cornerUpperLeftGreen[0] - _cornerUpperRightGreen[0] - _thicknessGreen );

    double dl    = 2 * ( horizontalLineWidth + verticalLineWidth ) / ( (double)_nProbingCellsLineGreen );
    double area  = dl * _thicknessGreen;
    double a_g   = 0;
    double l_max = _nProbingCellsLineGreen * dl;

    for( unsigned i = 0; i < _nProbingCellsLineGreen; i++ ) {    // Loop over probing cells
        _absorptionValsIntegrated[i] =
            ( _sigmaT[_probingCellsLineGreen[i]] - _sigmaS[_probingCellsLineGreen[i]] ) * scalarFlux[_probingCellsLineGreen[i]] * area;
        a_g += _absorptionValsIntegrated[i] / (double)_nProbingCellsLineGreen;
    }
    for( unsigned i = 0; i < _nProbingCellsLineGreen; i++ ) {    // Loop over probing cells
        _varAbsorptionValsIntegrated[i] = dl / l_max * ( a_g - _absorptionValsIntegrated[i] ) * ( a_g - _absorptionValsIntegrated[i] );
    }
}

std::vector<unsigned> SymmetricHohlraum::linspace2D( const std::vector<double>& start, const std::vector<double>& end, unsigned num_points ) {
    std::vector<unsigned> result;
    result.resize( num_points );
    double stepX = ( end[0] - start[0] ) / ( num_points - 1 );
    double stepY = ( end[1] - start[1] ) / ( num_points - 1 );

    for( unsigned i = 0; i < num_points; ++i ) {
        double x = start[0] + i * stepX;
        double y = start[1] + i * stepY;

        result[i] = _mesh->GetCellOfKoordinate( x, y );
    }

    return result;
}
// -------------- Moment Symmetric Hohlraum ---------------

SymmetricHohlraum_Moment::SymmetricHohlraum_Moment( Config* settings, Mesh* mesh, QuadratureBase* quad )
    : SymmetricHohlraum( settings, mesh, quad ) {}

SymmetricHohlraum_Moment::~SymmetricHohlraum_Moment() {}

std::vector<VectorVector> SymmetricHohlraum_Moment::GetExternalSource( const Vector& /* energies */ ) {
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

void SymmetricHohlraum_Moment::ComputeCurrentProbeMoment( const VectorVector& solution ) {
    for( unsigned idx_cell = 0; idx_cell < 4; idx_cell++ ) {    // Loop over probing cells
        _probingMoments[idx_cell][0] = solution[_probingCells[idx_cell]][0];
        if( _probingMoments[idx_cell].size() > 1 ) {
            _probingMoments[idx_cell][1] = solution[_probingCells[idx_cell]][1];
            _probingMoments[idx_cell][2] = solution[_probingCells[idx_cell]][2];
        }
    }
}