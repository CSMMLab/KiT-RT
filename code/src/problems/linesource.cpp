#include "problems/linesource.h"

// ---- LineSource_SN ----

LineSource_SN::LineSource_SN( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _physics = nullptr; }

LineSource_SN::~LineSource_SN(){};

VectorVector LineSource_SN::GetScatteringXS( const std::vector<double>& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 1.0 ) );
}

VectorVector LineSource_SN::GetTotalXS( const std::vector<double>& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 1.0 ) );
}

std::vector<VectorVector> LineSource_SN::GetExternalSource( const std::vector<double>& energies ) {
    return std::vector<VectorVector>( 1u, std::vector<Vector>( _mesh->GetNumCells(), Vector( 1u, 0.0 ) ) );
}

std::vector<double> LineSource_SN::GetStoppingPower( const std::vector<double>& energies ) {
    // @TODO
    return std::vector<double>( energies.size(), 0.0 );
}

VectorVector LineSource_SN::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids = _mesh->GetCellMidPoints();
    double t      = 3.2e-4;    // pseudo time for gaussian smoothing
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];
        psi[j]   = 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) );
    }
    return psi;
}

// ---- LineSource_PN ----

int LineSource_PN::GlobalIndex( int l, int k ) const {
    int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
    int prevIndicesThisLevel = k + l;    // number of previous indices in current level
    return numIndicesPrevLevel + prevIndicesThisLevel;
}

LineSource_PN::LineSource_PN( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _physics = nullptr; }

LineSource_PN::~LineSource_PN(){};

VectorVector LineSource_PN::GetScatteringXS( const std::vector<double>& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 1.0 ) );
}

VectorVector LineSource_PN::GetTotalXS( const std::vector<double>& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), 1.0 ) );
}

std::vector<VectorVector> LineSource_PN::GetExternalSource( const std::vector<double>& energies ) {
    return std::vector<VectorVector>( 1u, std::vector<Vector>( _mesh->GetNumCells(), Vector( 1u, 0.0 ) ) );
}

std::vector<double> LineSource_PN::GetStoppingPower( const std::vector<double>& energies ) {
    // @TODO
    return std::vector<double>( energies.size(), 0.0 );
}

VectorVector LineSource_PN::SetupIC() {
    // Compute number of equations in the system
    int ntotalEquations = GlobalIndex( _settings->GetMaxMomentDegree(), _settings->GetMaxMomentDegree() ) + 1;

    VectorVector psi( _mesh->GetNumCells(), Vector( ntotalEquations, 0 ) );    // zero could lead to problems?
    VectorVector cellMids = _mesh->GetCellMidPoints();

    // Initial condition is dirac impulse at (x,y) = (0,0) ==> constant in angle ==> all moments are zero.
    double t = 3.2e-4;    // pseudo time for gaussian smoothing (Approx to dirac impulse)
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x  = cellMids[j][0];
        double y  = cellMids[j][1];
        psi[j][0] = 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) );
    }

    // for( unsigned j = 0; j < cellMids.size(); ++j ) {
    //    if( cellMids[j][0] < 0 && cellMids[j][1] > 0 )
    //        psi[j][0] = 1.0;    // / ( 4.0 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) ) / ( 4 * M_PI );
    //    else
    //        psi[j][0] = 0.0;
    //}
    return psi;
}
