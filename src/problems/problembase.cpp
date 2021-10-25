#include "common/config.hpp"

#include "problems/aircavity1d.hpp"
#include "problems/checkerboard.hpp"
#include "problems/electronrt.hpp"
#include "problems/isotropicsource2d.hpp"
#include "problems/linesource.hpp"
#include "problems/musclebonelung.hpp"
#include "problems/phantom2d.hpp"
#include "problems/problembase.hpp"
#include "problems/waterphantom.hpp"

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
            if( settings->GetSolverName() == PN_SOLVER || settings->GetSolverName() == MN_SOLVER )
                return new LineSource_PN( settings, mesh );
            else
                return new LineSource_SN( settings, mesh );    // default
        }
        case PROBLEM_Checkerboard: {
            if( settings->GetSolverName() == PN_SOLVER || settings->GetSolverName() == MN_SOLVER )
                return new Checkerboard_Moment( settings, mesh );
            else
                return new Checkerboard_SN( settings, mesh );    // default
        }
        case PROBLEM_ElectronRT: return new ElectronRT( settings, mesh );
        case PROBLEM_AirCavity: return new AirCavity1D( settings, mesh );
        case PROBLEM_MuscleBoneLung: return new MuscleBoneLung( settings, mesh );
        case PROBLEM_WaterPhantom: return new WaterPhantom( settings, mesh );
        case PROBLEM_Phantom2D: return new Phantom2D( settings, mesh );
        case PROBLEM_LineSource_Pseudo_1D: return new LineSource_SN_Pseudo1D( settings, mesh );
        case PROBLEM_LineSource_Pseudo_1D_Physics: return new LineSource_SN_Pseudo1D_Physics( settings, mesh );
        case PROBLEM_IsotropicSource_2D: return new IsotropicSource2D( settings, mesh );
        default: return new ElectronRT( settings, mesh );    // Use RadioTherapy as dummy
    }
}

// Default densities = 1
std::vector<double> ProblemBase::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }

// Legacy code: Scattering crossection loaded from database ENDF with physics class -> later overwritten with ICRU data
VectorVector ProblemBase::GetScatteringXSE( const Vector& energies, const Vector& angles ) { return _physics->GetScatteringXS( energies, angles ); }

// Stopping powers from phyics class or default = -1
Vector ProblemBase::GetStoppingPower( const Vector& energies ) {
    if( _physics ) {
        return _physics->GetStoppingPower( energies );
    }
    else {
        ErrorMessages::Error( "Problem child class has not initialized a 'Physics' object!", CURRENT_FUNCTION );
        return Vector( 1, -1.0 );
    }
}
