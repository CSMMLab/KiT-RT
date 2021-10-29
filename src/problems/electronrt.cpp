#include "problems/electronrt.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"

// Constructor: Legacy code, physics class is no longer used
ElectronRT::ElectronRT( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {
    _physics = new EPICS( settings->GetHydrogenFile(), settings->GetOxygenFile(), settings->GetStoppingPowerFile() );
}

ElectronRT::~ElectronRT() { delete _physics; }

VectorVector ElectronRT::GetScatteringXS( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( 1, Vector( 1, 0.0 ) );
}

VectorVector ElectronRT::GetTotalXS( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( 1, Vector( 1, 0.0 ) );
}

std::vector<Matrix> ElectronRT::GetScatteringXSE( const Vector& /*energies*/, const Matrix& /*angles*/ ) {
    // @TODO
    // Specified in subclasses
    // return _physics->GetScatteringXS( energies, angles );
    return std::vector<Matrix>( 1, Matrix( 1, 1 ) );
}

Vector ElectronRT::GetTotalXSE( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    // return _physics->GetTotalXSE( energies );
    return Vector( 1 );
}

std::vector<VectorVector> ElectronRT::GetExternalSource( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    return std::vector<VectorVector>( 1, std::vector<Vector>( 1, Vector( 1, 0.0 ) ) );
}

VectorVector ElectronRT::SetupIC() {
    // @TODO
    // Specified in subclasses
    // return VectorVector( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    return VectorVector( 1, Vector( 1, 1e-10 ) );
}

void ElectronRT::LoadXSH20( std::string /*fileSigmaS*/, std::string /*fileSigmaT*/ ) {
    // @TODO
    // Specified in subclasses
}

// Default densities = 1, for patient files this is overwritten with "real" densities
std::vector<double> ElectronRT::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }
