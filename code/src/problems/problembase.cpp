#include "problems/problembase.h"
#include "problems/checkerboard.h"
#include "problems/electronrt.h"
#include "problems/linesource.h"

ProblemBase::ProblemBase( Config* settings, Mesh* mesh ) : _settings( settings ), _mesh( mesh ) {}

ProblemBase::~ProblemBase() {}

ProblemBase* ProblemBase::Create( Config* settings, Mesh* mesh ) {
    auto name = settings->GetProblemName();
    switch( name ) {
        case PROBLEM_LineSource: return new LineSource( settings, mesh );
        case PROBLEM_Checkerboard: return new Checkerboard( settings, mesh );
        case PROBLEM_ElectronRT: return new ElectronRT( settings, mesh );
        default: return new ElectronRT( settings, mesh );    // Use RadioTherapy as dummy
    }
}
