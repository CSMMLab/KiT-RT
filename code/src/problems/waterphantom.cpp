#include "problems/waterphantom.h"
#include "common/config.h"
#include "common/mesh.h"

WaterPhantom::WaterPhantom( Config* settings, Mesh* mesh ) : ElectronRT( settings, mesh ) {
    _physics = new Physics( settings->GetHydrogenFile(), settings->GetOxygenFile(), "../input/stopping_power.txt" );    // TODO
}

WaterPhantom::~WaterPhantom() { delete _physics; }

std::vector<VectorVector> WaterPhantom::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector WaterPhantom::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    double s      = 0.1;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        psi[j]   = 1.0 / ( s * sqrt( 2 * M_PI ) ) * std::exp( -x * x / ( 2 * s * s ) );
    }
    return psi;
}

std::vector<double> WaterPhantom::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }
