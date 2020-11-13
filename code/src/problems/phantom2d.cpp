#include "problems/phantom2d.h"
#include "common/config.h"
#include "common/mesh.h"

Phantom2D::Phantom2D( Config* settings, Mesh* mesh ) : ElectronRT( settings, mesh ) {}

Phantom2D::~Phantom2D() { delete _physics; }

std::vector<VectorVector> Phantom2D::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector Phantom2D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    double s      = 0.1;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x                                = cellMids[j][0];
        psi[j][_settings->GetNQuadPoints() - 1] = 1.0 / ( s * sqrt( 2 * M_PI ) ) * std::exp( -x * x / ( 2 * s * s ) );
    }
    return psi;
}

std::vector<double> Phantom2D::GetDensity( const VectorVector& cellMidPoints ) {
    std::string imageFile = "../tests/input/phantom.png";
    std::string meshFile  = _settings->GetMeshFile();
    Matrix gsImage = createSU2MeshFromImage( imageFile, meshFile );
    auto bounds      = _mesh->GetBounds();

    double xMin = bounds[0].first;
    double xMax = bounds[0].second;
    double yMin = bounds[1].first;
    double yMax = bounds[1].second;

    unsigned m = gsImage.rows();
    unsigned n = gsImage.columns();

    Vector x( m + 1 ), y( n + 1 );
    for( unsigned i = 0; i < m + 1; ++i ) {
        x[i] = static_cast<double>( i ) / static_cast<double>( m ) * ( xMax - xMin );
    }
    for( unsigned i = 0; i < n + 1; ++i ) y[i] = static_cast<double>( i ) / static_cast<double>( n ) * ( yMax - yMin );

    Interpolation interp( x, y, gsImage );
    std::vector<double> result( _mesh->GetNumCells(), 0.0 );
    for( unsigned i = 0; i < _mesh->GetNumCells(); ++i ) {
        result[i] = std::clamp( interp( cellMidPoints[i][0], cellMidPoints[i][1] ), 0.0, 1.0 );
    }

 
  }