#include "problems/checkerboard.h"

Checkerboard::Checkerboard( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {
    _physics = nullptr;
    _scatteringXS.resize( _mesh->GetNumCells(), 1.0 );
    _totalXS.resize( _mesh->GetNumCells(), 1.0 );
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isAbsorption( cellMids[j] ) ) {
            _scatteringXS[j] = 0.0;
            _totalXS[j]      = 10.0;
        }
    }
}

Checkerboard::~Checkerboard() {}

VectorVector Checkerboard::GetScatteringXS( const std::vector<double>& energies ) { return VectorVector( energies.size(), _scatteringXS ); }

VectorVector Checkerboard::GetTotalXS( const std::vector<double>& energies ) { return VectorVector( energies.size(), _totalXS ); }

VectorVector Checkerboard::GetExternalSource( const std::vector<double>& energies ) {
    std::vector<unsigned> marker;
    auto cellMids = _mesh->GetCellMidPoints();
    auto cellArea = _mesh->GetCellAreas();
    for( unsigned j = 0; j < _mesh->GetNumCells(); ++j ) {
        if( cellMids[j][0] >= 2.5 && cellMids[j][0] <= 3.5 && cellMids[j][1] >= 2.5 && cellMids[j][1] <= 3.5 ) {
            marker.push_back( j );
        }
    }
    VectorVector Q( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) );
    for( unsigned i = 0; i < energies.size(); ++i ) {
        for( auto j : marker ) {
            Q[i][j] = cellArea[j];
        }
    }
    return Q;
}

std::vector<double> Checkerboard::GetStoppingPower( const std::vector<double>& energies ) {
    // @TODO
    return std::vector<double>( energies.size(), 0.0 );
}

VectorVector Checkerboard::SetupIC() { return VectorVector( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-7 ) ); }

bool Checkerboard::isAbsorption( const Vector& pos ) {
    std::vector<double> lbounds{ 0.5, 1.5, 2.5, 3.5, 4.5 };
    std::vector<double> ubounds{ 1.5, 2.5, 3.5, 4.5, 5.5 };
    for( unsigned k = 0; k < lbounds.size(); ++k ) {
        for( unsigned l = 0; l < lbounds.size(); ++l ) {
            if( ( l + k ) % 2 == 1 || ( k == 2 && l == 2 ) || ( k == 4 && l == 2 ) ) continue;
            if( pos[0] >= lbounds[k] && pos[0] <= ubounds[k] && pos[1] >= lbounds[l] && pos[1] <= ubounds[l] ) {
                return true;
            }
        }
    }
    return false;
}
