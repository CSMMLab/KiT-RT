#include "common/config.h"

#include "problems/aircavity1d.h"
#include "problems/checkerboard.h"
#include "problems/electronrt.h"
#include "problems/linesource.h"
#include "problems/musclebonelung.h"
#include "problems/problembase.h"
#include "problems/waterphantom.h"

ProblemBase::ProblemBase( Config* settings, Mesh* mesh ) {
    _settings = settings;
    _mesh     = mesh;
}

ProblemBase::~ProblemBase() {}

ProblemBase* ProblemBase::Create( Config* settings, Mesh* mesh ) {
    auto name = settings->GetProblemName();
    switch( name ) {
        case PROBLEM_LineSource: {
            if( settings->GetSolverName() == PN_SOLVER || settings->GetSolverName() == MN_SOLVER )
                return new LineSource_PN( settings, mesh );
            else
                return new LineSource_SN( settings, mesh );    // default
        }
        case PROBLEM_Checkerboard: {
            if( settings->GetSolverName() == PN_SOLVER || settings->GetSolverName() == MN_SOLVER )
                return new Checkerboard_PN( settings, mesh );
            else
                return new Checkerboard_SN( settings, mesh );    // default
        }
        case PROBLEM_ElectronRT: return new ElectronRT( settings, mesh );
        case PROBLEM_AirCavity: return new AirCavity1D( settings, mesh );
        case PROBLEM_MuscleBoneLung: return new MuscleBoneLung( settings, mesh );
        case PROBLEM_WaterPhantom: return new WaterPhantom( settings, mesh );
        case PROBLEM_LineSource_Pseudo_1D: return new LineSource_SN_Pseudo1D( settings, mesh );
        case PROBLEM_LineSource_Pseudo_1D_Physics: return new LineSource_SN_Pseudo1D_Physics( settings, mesh );
        default: return new ElectronRT( settings, mesh );    // Use RadioTherapy as dummy
    }
}

std::vector<double> ProblemBase::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }

VectorVector ProblemBase::GetScatteringXSE( const Vector& energies, const Vector& angles ) { return _physics->GetScatteringXS( energies, angles ); }

Vector ProblemBase::GetStoppingPower( const Vector& energies ) {
    if( _physics ) {
        return _physics->GetStoppingPower( energies );
    }
    else {
        ErrorMessages::Error( "Problem child class has not initialized a 'Physics' object!", CURRENT_FUNCTION );
        return Vector( 1, -1.0 );
    }
}
