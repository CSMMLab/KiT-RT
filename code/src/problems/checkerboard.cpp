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
    Vector Q( _mesh->GetNumCells(), 0.0 );
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( cellMids[j][0] >= 3 && cellMids[j][0] <= 4 && cellMids[j][1] >= 3 && cellMids[j][1] <= 4 ) Q[j] = 1.0;
    }
    return VectorVector( energies.size(), Q );
}

std::vector<double> Checkerboard::GetStoppingPower( const std::vector<double>& energies ) {
    // @TODO
    return std::vector<double>( energies.size(), 0.0 );
}

VectorVector Checkerboard::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( cellMids[j][0] >= 3 && cellMids[j][0] <= 4 && cellMids[j][1] >= 3 && cellMids[j][1] <= 4 ) psi[j] = 1.0;
    }
    return psi;
}

bool Checkerboard::isAbsorption( const Vector& pos ) {
    std::vector<double> lbounds{ 1, 2, 3, 4, 5 };
    std::vector<double> ubounds{ 2, 3, 4, 5, 6 };
    for( unsigned k = 0; k < lbounds.size(); ++k ) {
        for( unsigned l = 0; l < lbounds.size(); ++l ) {
            if( ( l + k ) % 2 == 1 || ( k == 2 && l == 2 ) || ( k == 2 && l == 4 ) ) continue;
            if( pos[0] >= lbounds[k] && pos[0] <= ubounds[k] && pos[1] >= lbounds[l] && pos[1] <= ubounds[l] ) {
                return true;
            }
        }
    }
    return false;
}
