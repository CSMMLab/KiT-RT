#include "problems/aircavity1d.h"
#include "common/config.h"
#include "common/mesh.h"

AirCavity1D::AirCavity1D( Config* settings, Mesh* mesh ) : ElectronRT( settings, mesh ) {}

AirCavity1D::~AirCavity1D() { delete _physics; }

std::vector<VectorVector> AirCavity1D::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector AirCavity1D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    double s      = 0.1;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x                                = cellMids[j][0];
        psi[j][_settings->GetNQuadPoints() - 1] = 1.0 / ( s * sqrt( 2 * M_PI ) ) * std::exp( -x * x / ( 2 * s * s ) );
    }
    return psi;
}

std::vector<double> AirCavity1D::GetDensity( const VectorVector& cellMidPoints ) { 
    std::vector<double> densities ( 4*cellMidPoints.size()/9, 1.0 );
    std::vector<double> air( 2*cellMidPoints.size()/9, 0.00125  );
    std::vector<double> water ( 3*cellMidPoints.size()/9, 1.0 );
    densities.insert(densities.end(),air.begin(),air.end());
    densities.insert(densities.end(),water.begin(),water.end());
    return densities; }
