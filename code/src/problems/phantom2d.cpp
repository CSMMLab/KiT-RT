#include "problems/phantom2d.h"
#include "blaze/math/smp/default/DenseVector.h"    // for smpAssign
#include "common/config.h"
#include "common/io.h"    // for createSU2MeshFromImage
#include "common/mesh.h"
#include "problems/electronrt.h"        // for ElectronRT
#include "problems/epics.h"             // for EPICS
#include "toolboxes/interpolation.h"    // for Interpolation
#include <algorithm>                    // for clamp
#include <ext/alloc_traits.h>           // for __alloc_traits<>::va...
#include <math.h>                       // for exp, sqrt, M_PI
#include <memory>                       // for allocator_traits<>::...
#include <string>                       // for string

Phantom2D::Phantom2D( Config* settings, Mesh* mesh ) : ElectronRT( settings, mesh ) {}

Phantom2D::~Phantom2D() { delete _physics; }

std::vector<VectorVector> Phantom2D::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector Phantom2D::SetupIC() {
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

std::vector<double> Phantom2D::GetDensity( const VectorVector& cellMidPoints ) {

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
