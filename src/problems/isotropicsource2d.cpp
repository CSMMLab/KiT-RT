#include "problems/isotropicsource2d.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "problems/problembase.hpp"
#include "solvers/csdpn_starmap_constants.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/interpolation.hpp"
#include "toolboxes/sphericalbase.hpp"

#include <fstream>
#include <numeric>

IsotropicSource2D::IsotropicSource2D( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {}

IsotropicSource2D::~IsotropicSource2D() {}

VectorVector IsotropicSource2D::GetScatteringXS( const Vector& energies ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( 1, Vector( 1, 0.0 ) );
}

VectorVector IsotropicSource2D::GetTotalXS( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( 1, Vector( 1, 0.0 ) );
}

std::vector<Matrix> IsotropicSource2D::GetScatteringXSE( const Vector& /*energies*/, const Matrix& /*angles*/ ) {
    // @TODO
    // Specified in subclasses
    // return _physics->GetScatteringXS( energies, angles );
    return std::vector<Matrix>( 1, Matrix( 1, 1 ) );
}

Vector IsotropicSource2D::GetTotalXSE( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    // return _physics->GetTotalXSE( energies );
    return Vector( 1 );
}

std::vector<VectorVector> IsotropicSource2D::GetExternalSource( const Vector& energies ) {
    auto zeroVec = Vector( _settings->GetNQuadPoints(), 0.0 );
    auto uniform = std::vector<Vector>( _mesh->GetNumCells(), zeroVec );
    auto Q       = std::vector<VectorVector>( energies.size(), uniform );
    return Q;
}

VectorVector IsotropicSource2D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids         = _mesh->GetCellMidPoints();
    double enterPositionX = 0.5;    // 0.0;
    double enterPositionY = 0.5;
    // auto boundaryCells    = _mesh->GetBoundaryTypes();
    // Case 1: Ingoing radiation in just one cell
    // find cell that best matches enter position
    // double dist = 1000.0;
    // unsigned indexSource = 0;
    /*
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        // if( boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];
        if( x >= 0.49 && x <= 0.5 && y >= 0.49 && y <= 0.5 ) {
            psi[j] = Vector( _settings->GetNQuadPoints(), 1.0 );
        }
        //}
    }*/
    // psi[indexSource] = Vector( _settings->GetNQuadPoints(), 1.0 );

    // Case 2: Ingoing radiation as Gauss curve
    double t = 1e-5;    // pseudo time for gaussian smoothing
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0] - enterPositionX;
        double y = cellMids[j][1] - enterPositionY;
        psi[j]   = Vector( _settings->GetNQuadPoints(), 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) ) );
    }

    return psi;
}

std::vector<double> IsotropicSource2D::GetDensity( const VectorVector& /*cellMidPoints*/ ) {
    double rhoL = 1.0;
    double rhoR = 5.0;
    std::vector<double> rho( _settings->GetNCells(), rhoL );

    // use values rhoR on right third of domain
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        if( x >= 0.56 ) {
            rho[j] = rhoR;
        }
    }
    return rho;
}

// Moment  version below
IsotropicSource2D_Moment::IsotropicSource2D_Moment( Config* settings, Mesh* mesh ) : IsotropicSource2D( settings, mesh ) {}

IsotropicSource2D_Moment::~IsotropicSource2D_Moment() {}

VectorVector IsotropicSource2D_Moment::SetupIC() {

    SphericalBase* tempBase  = SphericalBase::Create( _settings );
    unsigned ntotalEquations = tempBase->GetBasisSize();
    delete tempBase;
    // write initial condition
    Vector pos_beam              = Vector{ 0.5, 0.5 };
    VectorVector initialSolution = VectorVector( _mesh->GetNumCells(), Vector( ntotalEquations, 0.0 ) );
    VectorVector cellMidpoints   = _mesh->GetCellMidPoints();
    const double stddev          = .01;
    double x = 0.0, y = 0.0, f = 0.0;

    for( unsigned idx_cell = 0; idx_cell < _mesh->GetNumCells(); ++idx_cell ) {
        x                            = cellMidpoints[idx_cell][0];
        y                            = cellMidpoints[idx_cell][1];
        f                            = NormPDF( x, pos_beam[0], stddev ) * NormPDF( y, pos_beam[1], stddev );
        initialSolution[idx_cell][0] = f * StarMAPmoments[0];
        for( unsigned idx_sys = 1; idx_sys < ntotalEquations; idx_sys++ ) {
            initialSolution[idx_cell][idx_sys] = f * StarMAPmoments[idx_sys];    // must be VectorVector
        }
    }
    return initialSolution;
}

double IsotropicSource2D_Moment::NormPDF( double x, double mu, double sigma ) {
    return INV_SQRT_2PI / sigma * std::exp( -( ( x - mu ) * ( x - mu ) ) / ( 2.0 * sigma * sigma ) );
}
