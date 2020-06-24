#include "problems/waterphantom.h"

WaterPhantom::WaterPhantom( Config* settings, Mesh* mesh ) : ElectronRT( settings, mesh ) {
    // @TODO get pointer to correct physics class
    _physics = nullptr;
}

WaterPhantom::~WaterPhantom() {}

std::vector<VectorVector> WaterPhantom::GetExternalSource( const std::vector<double>& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

std::vector<double> WaterPhantom::GetStoppingPower( const std::vector<double>& energies ) {
    // @TODO get correct stopping power
    return std::vector<double>( energies.size(), 1.0 );
}

VectorVector WaterPhantom::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    double t      = 3.2e-4;    // pseudo time for gaussian smoothing
    double s      = 0.1;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        psi[j]   = 1.0 / ( s * sqrt( 2 * M_PI ) ) * std::exp( -x * x / ( 2 * s * s ) );
    }
    return psi;
}

std::vector<double> WaterPhantom::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }
