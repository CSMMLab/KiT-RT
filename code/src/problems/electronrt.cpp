#include "problems/electronrt.h"

ElectronRT::ElectronRT( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {
    // @TODO
    _physics = nullptr;
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

VectorVector ElectronRT::GetExternalSource( const std::vector<double>& energies ) {
    // @TODO
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
}

std::vector<double> ElectronRT::GetStoppingPower( const std::vector<double>& energies ) {
    // @TODO
    return std::vector<double>( energies.size(), 0.0 );
}

VectorVector ElectronRT::SetupIC() {
    // @TODO
    return VectorVector( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
}

void ElectronRT::LoadXSH20( std::string fileSigmaS, std::string fileSigmaT ) {
    // @TODO
}
