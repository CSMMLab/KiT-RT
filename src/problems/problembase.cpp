#include "problems/problembase.hpp"

#include "common/config.hpp"
#include "common/mesh.hpp"
#include "problems/aircavity1d.hpp"
#include "problems/checkerboard.hpp"
#include "problems/hohlraum.hpp"
#include "problems/lattice.hpp"
#include "problems/linesource.hpp"
#include "problems/meltingcube.hpp"
#include "problems/phantomimage.hpp"
#include "problems/radiationctimage.hpp"
#include "problems/starmapvalidation.hpp"
#include "problems/symmetrichohlraum.hpp"
#include "toolboxes/errormessages.hpp"

ProblemBase::ProblemBase( Config* settings, Mesh* mesh ) {
    _settings = settings;
    _mesh     = mesh;
    SetGhostCells();
}

ProblemBase::~ProblemBase() {}

ProblemBase* ProblemBase::Create( Config* settings, Mesh* mesh ) {
    // Choose problem type
    switch( settings->GetProblemName() ) {
        case PROBLEM_Linesource: {
            if( settings->GetIsMomentSolver() )
                return new LineSource_Moment( settings, mesh );
            else
                return new LineSource_SN( settings, mesh );
        } break;
        case PROBLEM_Linesource1D: {
            if( settings->GetIsMomentSolver() )
                return new LineSource_Moment_1D( settings, mesh );
            else
                return new LineSource_SN_1D( settings, mesh );
        } break;
        case PROBLEM_Checkerboard: {
            if( settings->GetIsMomentSolver() )
                return new Checkerboard_Moment( settings, mesh );
            else
                return new Checkerboard_SN( settings, mesh );
        } break;
        case PROBLEM_Checkerboard1D: {
            if( settings->GetIsMomentSolver() )
                return new Checkerboard_Moment_1D( settings, mesh );
            else
                return new Checkerboard_SN_1D( settings, mesh );
        } break;
        case PROBLEM_Aircavity1D: {
            if( settings->GetIsMomentSolver() )
                return new AirCavity1D_Moment( settings, mesh );
            else
                return new AirCavity1D( settings, mesh );
        } break;
        case PROBLEM_StarmapValidation: {
            if( settings->GetIsMomentSolver() )
                return new StarMapValidation_Moment( settings, mesh );
            else
                return new StarMapValidation_SN( settings, mesh );
        } break;
        case PROBLEM_Phantomimage: return new PhantomImage( settings, mesh );
        case PROBLEM_RadiationCT: {
            if( settings->GetIsMomentSolver() )
                return new RadiationCTImage_Moment( settings, mesh );
            else
                return new RadiationCTImage( settings, mesh );
        } break;
        case PROBLEM_Meltingcube: {
            if( settings->GetIsMomentSolver() )
                return new MeltingCube_Moment( settings, mesh );
            else
                return new MeltingCube_SN( settings, mesh );
        } break;
        case PROBLEM_Meltingcube1D: {
            if( settings->GetIsMomentSolver() )
                return new MeltingCube_Moment_1D( settings, mesh );
            else
                return new MeltingCube_SN_1D( settings, mesh );
        } break;
        case PROBLEM_Hohlraum: {
            if( settings->GetIsMomentSolver() )
                return new Hohlraum_Moment( settings, mesh );
            else
                return new Hohlraum( settings, mesh );
        } break;
        case PROBLEM_SymmetricHohlraum: {
            if( settings->GetIsMomentSolver() )
                return new SymmetricHohlraum_Moment( settings, mesh );
            else
                return new SymmetricHohlraum( settings, mesh );
        } break;
        case PROBLEM_Lattice: {
            if( settings->GetIsMomentSolver() )
                return new Lattice_Moment( settings, mesh );
            else
                return new Lattice_SN( settings, mesh );
        } break;

        default: ErrorMessages::Error( "No valid physical problem chosen. Please check your config file", CURRENT_FUNCTION ); return nullptr;
    }
}

// Default densities = 1
std::vector<double> ProblemBase::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }

// Legacy code: Scattering crossection loaded from database ENDF with physics
// class -> later overwritten with ICRU data
VectorVector ProblemBase::GetScatteringXSE( const Vector& /*energies*/, const Vector& /*angles*/ ) {
    ErrorMessages::Error( "Not yet implemented", CURRENT_FUNCTION );
    return VectorVector( 1, Vector( 1, 0 ) );
}

// Stopping powers from phyics class or default = -1
Vector ProblemBase::GetStoppingPower( const Vector& /* energies */ ) {
    ErrorMessages::Error( "Not yet implemented", CURRENT_FUNCTION );
    return Vector( 1, -1.0 );
}

void ProblemBase::SetGhostCells() {
    // Loop over all cells. If its a Dirichlet boundary, add cell to dict with
    // {cell_idx, boundary_value}
    auto cellBoundaries = _mesh->GetBoundaryTypes();
    std::map<int, Vector> ghostCellMap;

    Vector dummyGhostCell( 1, 0.0 );

    for( unsigned idx_cell = 0; idx_cell < _mesh->GetNumCells(); idx_cell++ ) {
        if( cellBoundaries[idx_cell] == BOUNDARY_TYPE::NEUMANN || cellBoundaries[idx_cell] == BOUNDARY_TYPE::DIRICHLET ) {
            // TODO: Refactor Boundary Conditions: We only have Ghost Cells with
            // Dirichlet conditions right now
            ghostCellMap.insert( { idx_cell, dummyGhostCell } );
        }
        _ghostCells = ghostCellMap;
    }
}

const Vector& ProblemBase::GetGhostCellValue( int idx_cell, const Vector& cell_sol ) { return cell_sol; }
