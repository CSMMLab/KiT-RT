#include "problems/electronrt.h"
#include "common/config.h"
#include "common/mesh.h"

ElectronRT::ElectronRT( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {
    // @TODO
    _physics = new Physics();
}

ElectronRT::~ElectronRT() {}

VectorVector ElectronRT::GetScatteringXS( const std::vector<double>& energies ) {
    // @TODO
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
}

VectorVector ElectronRT::GetTotalXS( const std::vector<double>& energies ) {
    // @TODO
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
}

VectorVector ElectronRT::GetScatteringXSE( const Vector& energies, const Vector& angles ) {
    // @TODO
    _physics->GetScatteringXS( energies, angles );
    return VectorVector( energies.size(), Vector( angles.size(), 0.0 ) );
}

Vector ElectronRT::GetTotalXSE( const Vector& energies ) {
    // @TODO
    return Vector( energies.size(), 0.0 );
}

std::vector<VectorVector> ElectronRT::GetExternalSource( const std::vector<double>& energies ) {
    // @TODO
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

std::vector<double> ElectronRT::GetStoppingPower( const std::vector<double>& energies ) {
    // @TODO
    return std::vector<double>( energies.size(), 1.0 );
}

VectorVector ElectronRT::SetupIC() {
    // @TODO
    return VectorVector( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
}

void ElectronRT::LoadXSH20( std::string fileSigmaS, std::string fileSigmaT ) {
    // @TODO
}

std::vector<double> ElectronRT::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }
