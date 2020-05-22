#include "problems/electronrt.h"

ElectronRT::ElectronRT( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {
    // @TODO
    _physics = nullptr;
}

ElectronRT::~ElectronRT() {}

std::vector<Matrix> ElectronRT::GetScatteringXS( const double energy ) {
    // @TODO
    return std::vector<Matrix>( _mesh->GetNumCells(), Matrix( _settings->GetNQuadPoints(), _settings->GetNQuadPoints(), 0.0 ) );
}

std::vector<double> ElectronRT::GetTotalXS( const double energy ) {
    // @TODO
    return std::vector<double>( _mesh->GetNumCells(), 0.0 );
}

std::vector<double> ElectronRT::GetStoppingPower( const std::vector<double>& energies ) {
    // @TODO
    return std::vector<double>( energies.size(), 0.0 );
}

VectorVector ElectronRT::SetupIC() {
    // @TODO
    return VectorVector( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-7 ) );
}

void ElectronRT::LoadXSH20( std::string fileSigmaS, std::string fileSigmaT ) {
    // @TODO
}
