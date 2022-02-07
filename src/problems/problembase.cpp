#include "problems/problembase.hpp"
#include "common/config.hpp"
#include "problems/aircavity1d.hpp"
#include "problems/checkerboard.hpp"
#include "problems/linesource.hpp"
#include "problems/starmapvalidation.hpp"
#include "problems/waterphantom.hpp"
#include "toolboxes/errormessages.hpp"

ProblemBase::ProblemBase( Config* settings, Mesh* mesh ) {
    _settings = settings;
    _mesh     = mesh;
}

ProblemBase::~ProblemBase() {}

ProblemBase* ProblemBase::Create( Config* settings, Mesh* mesh ) {
    auto name = settings->GetProblemName();

    // Choose problem type
    switch( name ) {
        case PROBLEM_LineSource: {
            if( settings->GetSolverName() == PN_SOLVER || settings->GetSolverName() == MN_SOLVER ||
                settings->GetSolverName() == MN_SOLVER_NORMALIZED )
                return new LineSource_PN( settings, mesh );
            else
                return new LineSource_SN( settings, mesh );    // default
        }
        case PROBLEM_Checkerboard: {
            if( settings->GetSolverName() == PN_SOLVER || settings->GetSolverName() == MN_SOLVER ||
                settings->GetSolverName() == MN_SOLVER_NORMALIZED )
                return new Checkerboard_Moment( settings, mesh );
            else
                return new Checkerboard_SN( settings, mesh );    // default
        }
        case PROBLEM_AirCavity:
            if( settings->GetSolverName() == PN_SOLVER || settings->GetSolverName() == MN_SOLVER ||
                settings->GetSolverName() == MN_SOLVER_NORMALIZED )
                return new AirCavity1D_Moment( settings, mesh );
            else
                return new AirCavity1D( settings, mesh );    // default
        case PROBLEM_WaterPhantom: return new WaterPhantom1D( settings, mesh );
        case PROBLEM_Phantom2D: return new WaterPhantom( settings, mesh );
        case PROBLEM_LineSource_Pseudo_1D: return new LineSource_SN_Pseudo1D( settings, mesh );
        case PROBLEM_LinesourceDualDenstiy:
            if( settings->GetSolverName() == CSD_PN_SOLVER || settings->GetSolverName() == CSD_MN_SOLVER )
                return new StarMapValidation_Moment( settings, mesh );
            else                                                      // CSD_SN_SOLVER
                return new StarMapValidation_SN( settings, mesh );    // default
        default: ErrorMessages::Error( "No valid physical problem chosen. Please check your config file", CURRENT_FUNCTION ); return nullptr;
    }
}

// Default densities = 1
std::vector<double> ProblemBase::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }

// Legacy code: Scattering crossection loaded from database ENDF with physics class -> later overwritten with ICRU data
VectorVector ProblemBase::GetScatteringXSE( const Vector& energies, const Vector& angles ) {
    ErrorMessages::Error( "Not yet implemented", CURRENT_FUNCTION );
}

// Stopping powers from phyics class or default = -1
Vector ProblemBase::GetStoppingPower( const Vector& energies ) {
    ErrorMessages::Error( "Not yet implemented", CURRENT_FUNCTION );
    return Vector( 1, -1.0 );
}
