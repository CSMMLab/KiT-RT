#include "problems/electronrt.h"
#include "common/config.h"
#include "common/mesh.h"

// Constructor: Legacy code, physics class is no longer used
ElectronRT::ElectronRT( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {
    _physics = new Physics( settings->GetHydrogenFile(), settings->GetOxygenFile(), settings->GetStoppingPowerFile() );
}

ElectronRT::~ElectronRT() { delete _physics; }

VectorVector ElectronRT::GetScatteringXS( const Vector& energies ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
}

VectorVector ElectronRT::GetTotalXS( const Vector& energies ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
}

std::vector<Matrix> ElectronRT::GetScatteringXSE( const Vector& energies, const Matrix& angles ) {
    // @TODO
    // Specified in subclasses
    return _physics->GetScatteringXS( energies, angles );
}

Vector ElectronRT::GetTotalXSE( const Vector& energies ) {
    // @TODO
    // Specified in subclasses
    return _physics->GetTotalXSE( energies );
}

std::vector<VectorVector> ElectronRT::GetExternalSource( const Vector& energies ) {
    // @TODO
    // Specified in subclasses
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector ElectronRT::SetupIC() {
    // @TODO
    // Specified in subclasses
    return VectorVector( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
}

void ElectronRT::LoadXSH20( std::string /*fileSigmaS*/, std::string /*fileSigmaT*/ ) {
    // @TODO
    // Specified in subclasses
}

// Default densities = 1, for patient files this is overwritten with "real" densities
std::vector<double> ElectronRT::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }
