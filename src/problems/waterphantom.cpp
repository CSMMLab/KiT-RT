#include "problems/waterphantom.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "toolboxes/interpolation.hpp"

WaterPhantom1D::WaterPhantom1D( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _sigmaS = settings->GetSigmaS(); }

WaterPhantom1D::~WaterPhantom1D() {}

std::vector<VectorVector> WaterPhantom1D::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector WaterPhantom1D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) );
    return psi;
}

// Density of water = 1 everywhere
std::vector<double> WaterPhantom1D::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }

VectorVector WaterPhantom1D::GetScatteringXS( const Vector& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) );
}

VectorVector WaterPhantom1D::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) ); }

// ---- 2d test case

WaterPhantom::WaterPhantom( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {}

WaterPhantom::~WaterPhantom() {}

std::vector<VectorVector> WaterPhantom::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector WaterPhantom::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids         = _mesh->GetCellMidPoints();
    double s              = 0.1;
    double enterPositionX = 0.0;
    double enterPositionY = 0.5;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0] - enterPositionX;
        double y = cellMids[j][1] - enterPositionY;
        psi[j]   = 1.0 / ( s * sqrt( 2 * M_PI ) ) * std::exp( -( x * x + y * y ) / ( 2 * s * s ) );
    }
    return psi;
}

std::vector<double> WaterPhantom::GetDensity( const VectorVector& cellMidPoints ) {

    std::string imageFile = _settings->GetCTFile();
    std::string meshFile  = _settings->GetMeshFile();
    Matrix gsImage        = createSU2MeshFromImage( imageFile, meshFile );
    auto bounds           = _mesh->GetBounds();

    double xMin = bounds[0].first;
    double xMax = bounds[0].second;
    double yMin = bounds[1].first;
    double yMax = bounds[1].second;

    unsigned m = gsImage.rows();
    unsigned n = gsImage.columns();

    Vector x( m ), y( n );
    for( unsigned i = 0; i < m; ++i ) {
        x[i] = xMin + static_cast<double>( i ) / static_cast<double>( m - 1 ) * ( xMax - xMin );
    }
    for( unsigned i = 0; i < n; ++i ) y[i] = yMin + static_cast<double>( i ) / static_cast<double>( n - 1 ) * ( yMax - yMin );

    Interpolation interp( x, y, gsImage );
    std::vector<double> result( _mesh->GetNumCells(), 0.0 );
    for( unsigned i = 0; i < _mesh->GetNumCells(); ++i ) {
        result[i] = std::clamp( interp( cellMidPoints[i][0], cellMidPoints[i][1] ), 0.0, 1.85 );
    }
    return result;
}

VectorVector WaterPhantom::GetScatteringXS( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( 1, Vector( 1, 0.0 ) );
}

VectorVector WaterPhantom::GetTotalXS( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( 1, Vector( 1, 0.0 ) );
}
