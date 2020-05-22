#include "problems/linesource.h"

LineSource::LineSource( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _physics = nullptr; }

LineSource::~LineSource(){};

VectorVector LineSource::GetScatteringXS( const std::vector<double>& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
}

VectorVector LineSource::GetTotalXS( const std::vector<double>& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
}

std::vector<double> LineSource::GetStoppingPower( const std::vector<double>& energies ) {
    // @TODO
    return std::vector<double>( energies.size(), 0.0 );
}

VectorVector LineSource::SetupIC() {
    VectorVector psi = std::vector( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-7 ) );
    auto cellMids    = _mesh->GetCellMidPoints();
    Vector midPoint( 2, 0.5 );
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( norm( cellMids[j] - midPoint ) <= 0.1 ) {
            psi[j] = 1.0;
        }
    }
    return psi;
}
