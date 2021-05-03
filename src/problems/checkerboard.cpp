#include "problems/checkerboard.h"
#include "common/config.h"
#include "common/mesh.h"
#include "toolboxes/sphericalbase.h"

// ---- Checkerboard Sn ----
// Constructor for Ckeckerboard case with Sn
Checkerboard_SN::Checkerboard_SN( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {
    _physics = nullptr;

    // Initialise crosssections to 1
    _scatteringXS = Vector( _mesh->GetNumCells(), 1.0 );
    _totalXS      = Vector( _mesh->GetNumCells(), 1.0 );

    // For absorption cells: set scattering XS to 0 and absorption to 10
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isAbsorption( cellMids[j] ) ) {
            _scatteringXS[j] = 0.0;
            _totalXS[j]      = 10.0;
        }
    }
}

Checkerboard_SN::~Checkerboard_SN() {}

VectorVector Checkerboard_SN::GetScatteringXS( const Vector& energies ) { return VectorVector( energies.size(), _scatteringXS ); }

VectorVector Checkerboard_SN::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), _totalXS ); }

std::vector<VectorVector> Checkerboard_SN::GetExternalSource( const Vector& /*energies*/ ) {
    VectorVector Q( _mesh->GetNumCells(), Vector( 1u, 0.0 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isSource( cellMids[j] ) ) Q[j] = 1.0 / ( 4 * M_PI );    // isotropic source
    }
    return std::vector<VectorVector>( 1u, Q );
}

VectorVector Checkerboard_SN::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    return psi;
}

bool Checkerboard_SN::isAbsorption( const Vector& pos ) const {
    // Check whether pos is inside absorbing squares
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

bool Checkerboard_SN::isSource( const Vector& pos ) const {
    // Check whether pos is part of source region
    if( pos[0] >= 3 && pos[0] <= 4 && pos[1] >= 3 && pos[1] <= 4 )
        return true;
    else
        return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// ---- Checkerboard Moments ----

// Constructor for checkerboard case with Pn
Checkerboard_Moment::Checkerboard_Moment( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {

    _basis = SphericalBase::Create(_settings);
    _physics = nullptr;

    // Initialise crosssections = 1 (scattering)
    _scatteringXS = Vector( _mesh->GetNumCells(), 1.0 );
    _totalXS      = Vector( _mesh->GetNumCells(), 1.0 );

    // for absorption regions change crosssections to all absorption
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isAbsorption( cellMids[j] ) ) {
            _scatteringXS[j] = 0.0;
            _totalXS[j]      = 10.0;
        }
    }
}

Checkerboard_Moment::~Checkerboard_Moment() {
    delete _basis;
}

VectorVector Checkerboard_Moment::GetScatteringXS( const Vector& energies ) { return VectorVector( energies.size(), _scatteringXS ); }

VectorVector Checkerboard_Moment::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), _totalXS ); }

std::vector<VectorVector> Checkerboard_Moment::GetExternalSource( const Vector& /*energies*/ ) {
    VectorVector Q( _mesh->GetNumCells(), Vector( 1u, 0.0 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        if( isSource( cellMids[j] ) ) Q[j] = 1.0 / std::sqrt( 4 * M_PI );    // isotropic source
    }
    return std::vector<VectorVector>( 1u, Q );
}

VectorVector Checkerboard_Moment::SetupIC() {
    int ntotalEquations = _basis->GetBasisSize();
    VectorVector psi( _mesh->GetNumCells(), Vector( ntotalEquations, 1e-10 ) );
    return psi;
}

bool Checkerboard_Moment::isAbsorption( const Vector& pos ) const {
    // Check whether pos is in absorption region
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

bool Checkerboard_Moment::isSource( const Vector& pos ) const {
    // Check whether pos is in source region
    if( pos[0] >= 3 && pos[0] <= 4 && pos[1] >= 3 && pos[1] <= 4 )
        return true;
    else
        return false;
}
