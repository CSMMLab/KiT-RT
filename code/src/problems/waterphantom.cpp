#include "problems/waterphantom.h"
#include "blaze/math/smp/default/DenseVector.h"    // for smpAssign
#include "common/config.h"
#include "common/mesh.h"
#include "problems/electronrt.h"    // for ElectronRT

WaterPhantom::WaterPhantom( Config* settings, Mesh* mesh ) : ElectronRT( settings, mesh ) {}

WaterPhantom::~WaterPhantom() {}

std::vector<VectorVector> WaterPhantom::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector WaterPhantom::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) );
    return psi;
}

// Density of water = 1 everywhere
std::vector<double> WaterPhantom::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }
