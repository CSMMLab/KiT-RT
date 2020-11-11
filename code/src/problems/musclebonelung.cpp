#include "problems/musclebonelung.h"
#include "common/config.h"
#include "common/mesh.h"

MuscleBoneLung::MuscleBoneLung( Config* settings, Mesh* mesh ) : ElectronRT( settings, mesh ) {}

MuscleBoneLung::~MuscleBoneLung() { delete _physics; }

std::vector<VectorVector> MuscleBoneLung::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector MuscleBoneLung::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    double s      = 0.1;
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x                                = cellMids[j][0];
        psi[j][_settings->GetNQuadPoints() - 1] = 1.0 / ( s * sqrt( 2 * M_PI ) ) * std::exp( -x * x / ( 2 * s * s ) );
    }
    return psi;
}

std::vector<double> MuscleBoneLung::GetDensity( const VectorVector& cellMidPoints ) { 
    std::vector<double> densities ( 1.5*cellMidPoints.size()/9, 1.05 ); //muscle layer
    std::vector<double> bone( 1.5*cellMidPoints.size()/9, 1.92  );
    std::vector<double> lung ( 6*cellMidPoints.size()/9, 0.26 );
    densities.insert(densities.end(),bone.begin(),bone.end()-1);
    densities.insert(densities.end(),lung.begin(),lung.end()-1);
    return densities; }
