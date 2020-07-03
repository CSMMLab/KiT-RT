#include "problems/problembase.h"
#include "problems/checkerboard.h"
#include "problems/electronrt.h"
#include "problems/linesource.h"
#include "problems/waterphantom.h"

ProblemBase::ProblemBase( Config* settings, Mesh* mesh ) : _settings( settings ), _mesh( mesh ) {}

ProblemBase::~ProblemBase() {}

ProblemBase* ProblemBase::Create( Config* settings, Mesh* mesh ) {
    auto name = settings->GetProblemName();
    switch( name ) {
        case PROBLEM_LineSource: {
            if( settings->GetSolverName() == PN_SOLVER )
                return new LineSource_PN( settings, mesh );
            else
                return new LineSource_SN( settings, mesh );    // default
        }
        case PROBLEM_Checkerboard: return new Checkerboard( settings, mesh );
        case PROBLEM_ElectronRT: return new ElectronRT( settings, mesh );
        case PROBLEM_WaterPhantom: return new WaterPhantom( settings, mesh );
        case PROBLEM_LineSource_Pseudo_1D: return new LineSource_SN_Pseudo1D( settings, mesh );
        default: return new ElectronRT( settings, mesh );    // Use RadioTherapy as dummy
    }
}

std::vector<double> ProblemBase::GetDensity( const VectorVector& cellMidPoints ) { return std::vector<double>( cellMidPoints.size(), 1.0 ); }
