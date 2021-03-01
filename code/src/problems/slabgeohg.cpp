#include "problems/slabgeohg.h"
#include "blaze/math/smp/default/DenseVector.h"    // for smpAssign
#include "common/config.h"
#include "common/mesh.h"
#include "problems/problembase.h"    // for ProblemBase

// ---- SlabGeoHG ----

SlabGeoHG::SlabGeoHG( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _physics = nullptr; }

SlabGeoHG::~SlabGeoHG() {}

VectorVector SlabGeoHG::GetScatteringXS( const Vector& energies ) { return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) ); }

VectorVector SlabGeoHG::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 0.0 ) ); }

std::vector<VectorVector> SlabGeoHG::GetExternalSource( const Vector& /*energies*/ ) {
    return std::vector<VectorVector>( 1u, std::vector<Vector>( _mesh->GetNumCells(), Vector( 1u, 0.0 ) ) );
}

VectorVector SlabGeoHG::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) );
    // auto boundaryCells = _mesh->GetBoundaryTypes();
    // auto cellMids      = _mesh->GetCellMidPoints();
    // double t           = 3.2e-4;    // pseudo time for gaussian smoothing

    return psi;
}
