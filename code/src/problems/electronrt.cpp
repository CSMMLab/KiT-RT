#include "problems/electronrt.h"
#include "common/config.h"
#include "common/mesh.h"

ElectronRT::ElectronRT( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {
    _physics = new Physics( settings->GetHydrogenFile(), settings->GetOxygenFile(), "../input/stopping_power.txt" );    // TODO
}

ElectronRT::~ElectronRT() {}

VectorVector ElectronRT::GetScatteringXS( const Vector& energies ) {
    // @TODO
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
}

VectorVector ElectronRT::GetTotalXS( const Vector& energies ) {
    // @TODO
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
}

VectorVector ElectronRT::GetScatteringXSE( const Vector& energies, const Vector& angles ) {
    // @TODO
    std::cout << "before interpolation" << std::endl;
    _physics->GetScatteringXS( energies, angles );
    std::cout << "after interpolation" << std::endl;
    return VectorVector( energies.size(), Vector( angles.size(), 0.0 ) );
}

Vector ElectronRT::GetTotalXSE( const Vector& energies ) {
    // @TODO
    return Vector( energies.size(), 0.0 );
}

std::vector<VectorVector> ElectronRT::GetExternalSource( const Vector& energies ) {
    // @TODO
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

Vector ElectronRT::GetStoppingPower( const Vector& energies ) {
    // @TODO
    return Vector( energies.size(), 1.0 );
}

VectorVector ElectronRT::SetupIC() {
    // @TODO
    return VectorVector( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
}

void ElectronRT::LoadXSH20( std::string fileSigmaS, std::string fileSigmaT ) {
    // @TODO
}

std::vector<double> ElectronRT::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }
