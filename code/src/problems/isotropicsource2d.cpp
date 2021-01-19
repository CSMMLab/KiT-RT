#include "problems/isotropicsource2d.h"
#include "common/config.h"
#include "common/mesh.h"

IsotropicSource2D::IsotropicSource2D( Config* settings, Mesh* mesh ) : ElectronRT( settings, mesh ) {}

IsotropicSource2D::~IsotropicSource2D() { delete _physics; }

std::vector<VectorVector> IsotropicSource2D::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector IsotropicSource2D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids         = _mesh->GetCellMidPoints();
    double enterPositionX = 0.0;
    double enterPositionY = 0.5;
    auto boundaryCells    = _mesh->GetBoundaryTypes();

    // find cell that best matches enter position
    double dist          = 1000.0;
    unsigned indexSource = 0;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) {
            double x = cellMids[j][0] - enterPositionX;
            double y = cellMids[j][1] - enterPositionY;
            if( sqrt( x * x + y * y ) < dist ) {
                dist        = sqrt( x * x + y * y );
                indexSource = j;
            }
        }
    }
    psi[indexSource] = Vector( _settings->GetNQuadPoints(), 1.0 );
    return psi;
}

std::vector<double> IsotropicSource2D::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( _settings->GetNCells(), 1.0 ); }
