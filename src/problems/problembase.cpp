#include "problems/problembase.hpp"
#include "common/config.hpp"
#include "problems/aircavity1d.hpp"
#include "problems/checkerboard.hpp"
#include "problems/linesource.hpp"
#include "problems/phantomimage.hpp"
#include "problems/radiationctimage.hpp"
#include "problems/starmapvalidation.hpp"
#include "toolboxes/errormessages.hpp"

ProblemBase::ProblemBase( Config* settings, Mesh* mesh ) {
    _settings = settings;
    _mesh     = mesh;
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
