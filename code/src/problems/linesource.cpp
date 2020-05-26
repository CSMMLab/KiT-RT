#include "problems/linesource.h"

LineSource::LineSource( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _physics = nullptr; }

LineSource::~LineSource(){};

VectorVector LineSource::GetScatteringXS( const std::vector<double>& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 1.0 ) );
}

VectorVector LineSource::GetTotalXS( const std::vector<double>& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 1.0 ) );
}

VectorVector LineSource::GetExternalSource( const std::vector<double>& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
}

std::vector<double> LineSource::GetStoppingPower( const std::vector<double>& energies ) {
    // @TODO
    return std::vector<double>( energies.size(), 0.0 );
}

VectorVector LineSource::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    double t      = 3.2e-4;    // pseudo time for gaussian smoothing
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];
        psi[j]   = 1 / ( 4 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) ) / ( 4 * M_PI );
    }
    return psi;
}
