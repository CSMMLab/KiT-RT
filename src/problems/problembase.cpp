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
#include "problems/quarterhohlraum.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/errormessages.hpp"

ProblemBase::ProblemBase( Config* settings, Mesh* mesh, QuadratureBase* quad ) {
    _settings = settings;
    _mesh     = mesh;
    _quad     = quad;

    // initialize QOI helper variables
    _curMaxOrdinateOutflow   = 0.0;
    _curScalarOutflow        = 0.0;
    _totalScalarOutflow      = 0.0;
    _mass                    = 0.0;
    _changeRateFlux          = 0.0;
    _dummyProbeMoments       = VectorVector( 4, Vector( 3, 0.0 ) );
    _dummyProbeValsGreenLine = std::vector<double>( _settings->GetNumProbingCellsLineHohlraum(), 0.0 );
    SetGhostCells();
}

ProblemBase::~ProblemBase() {}

ProblemBase* ProblemBase::Create( Config* settings, Mesh* mesh, QuadratureBase* quad ) {
    // Choose problem type
    switch( settings->GetProblemName() ) {
        case PROBLEM_Linesource: {
            if( settings->GetIsMomentSolver() )
                return new LineSource_Moment( settings, mesh, quad );
            else
                return new LineSource_SN( settings, mesh, quad );
        } break;
        case PROBLEM_Linesource1D: {
            if( settings->GetIsMomentSolver() )
                return new LineSource_Moment_1D( settings, mesh, quad );
            else
                return new LineSource_SN_1D( settings, mesh, quad );
        } break;
        case PROBLEM_Checkerboard: {
            if( settings->GetIsMomentSolver() )
                return new Checkerboard_Moment( settings, mesh, quad );
            else
                return new Checkerboard_SN( settings, mesh, quad );
        } break;
        case PROBLEM_Checkerboard1D: {
            if( settings->GetIsMomentSolver() )
                return new Checkerboard_Moment_1D( settings, mesh, quad );
            else
                return new Checkerboard_SN_1D( settings, mesh, quad );
        } break;
        case PROBLEM_Aircavity1D: {
            if( settings->GetIsMomentSolver() )
                return new AirCavity1D_Moment( settings, mesh, quad );
            else
                return new AirCavity1D( settings, mesh, quad );
        } break;
        case PROBLEM_StarmapValidation: {
            if( settings->GetIsMomentSolver() )
                return new StarMapValidation_Moment( settings, mesh, quad );
            else
                return new StarMapValidation_SN( settings, mesh, quad );
        } break;
        case PROBLEM_Phantomimage: return new PhantomImage( settings, mesh, quad );
        case PROBLEM_RadiationCT: {
            if( settings->GetIsMomentSolver() )
                return new RadiationCTImage_Moment( settings, mesh, quad );
            else
                return new RadiationCTImage( settings, mesh, quad );
        } break;
        case PROBLEM_Meltingcube: {
            if( settings->GetIsMomentSolver() )
                return new MeltingCube_Moment( settings, mesh, quad );
            else
                return new MeltingCube_SN( settings, mesh, quad );
        } break;
        case PROBLEM_Meltingcube1D: {
            if( settings->GetIsMomentSolver() )
                return new MeltingCube_Moment_1D( settings, mesh, quad );
            else
                return new MeltingCube_SN_1D( settings, mesh, quad );
        } break;
        case PROBLEM_Hohlraum: {
            if( settings->GetIsMomentSolver() )
                return new Hohlraum_Moment( settings, mesh, quad );
            else
                return new Hohlraum( settings, mesh, quad );
        } break;
        case PROBLEM_SymmetricHohlraum: {
            if( settings->GetIsMomentSolver() )
                return new SymmetricHohlraum_Moment( settings, mesh, quad );
            else
                return new SymmetricHohlraum( settings, mesh, quad );
        } break;
        case PROBLEM_QuarterHohlraum: {
            if( settings->GetIsMomentSolver() )
                return new QuarterHohlraum_Moment( settings, mesh, quad );
            else
                return new QuarterHohlraum( settings, mesh, quad );
            
        } break;
        case PROBLEM_Lattice: {
            if( settings->GetIsMomentSolver() )
                return new Lattice_Moment( settings, mesh, quad );
            else
                return new Lattice_SN( settings, mesh, quad );
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

const Vector& ProblemBase::GetGhostCellValue( int /*idx_cell*/, const Vector& cell_sol ) { return cell_sol; }

void ProblemBase::ComputeMaxOrdinatewiseOutflow( const VectorVector& solution ) {
    if( _settings->GetSolverName() == SN_SOLVER || _settings->GetSolverName() == CSD_SN_SOLVER ) {
        double currOrdinatewiseOutflow = 0.0;
        double transportDirection      = 0.0;

        auto nCells   = _mesh->GetNumCells();
        auto cellMids = _mesh->GetCellMidPoints();
        auto areas    = _mesh->GetCellAreas();
        auto neigbors = _mesh->GetNeighbours();
        auto normals  = _mesh->GetNormals();

        auto quadPoints = _quad->GetPoints();
        auto weights    = _quad->GetWeights();
        auto nq         = _quad->GetNq();

        // Iterate over boundaries
        for( std::map<int, Vector>::iterator it = _ghostCells.begin(); it != _ghostCells.end(); ++it ) {
            int idx_cell = it->first;    // Get Boundary cell index
            for( unsigned idx_nbr = 0; idx_nbr < neigbors[idx_cell].size(); ++idx_nbr ) {
                // Find face that points outward
                if( neigbors[idx_cell][idx_nbr] == nCells ) {
                    // Iterate over transport directions
                    for( unsigned idx_quad = 0; idx_quad < nq; ++idx_quad ) {
                        transportDirection =
                            normals[idx_cell][idx_nbr][0] * quadPoints[idx_quad][0] + normals[idx_cell][idx_nbr][1] * quadPoints[idx_quad][1];
                        // Find outward facing transport directions
                        if( transportDirection > 0.0 ) {

                            currOrdinatewiseOutflow = transportDirection / norm( normals[idx_cell][idx_nbr] ) * solution[idx_cell][idx_quad];

                            if( currOrdinatewiseOutflow > _curMaxOrdinateOutflow ) _curMaxOrdinateOutflow = currOrdinatewiseOutflow;
                        }
                    }
                }
            }
        }
    }
    // TODO define alternative for Moment solvers
}

void ProblemBase::ComputeCurrentOutflow( const VectorVector& solution ) {
    if( _settings->GetSolverName() == SN_SOLVER || _settings->GetSolverName() == CSD_SN_SOLVER ) {

        _curScalarOutflow         = 0.0;
        double transportDirection = 0.0;

        auto nCells   = _mesh->GetNumCells();
        auto cellMids = _mesh->GetCellMidPoints();
        auto areas    = _mesh->GetCellAreas();
        auto neigbors = _mesh->GetNeighbours();
        auto normals  = _mesh->GetNormals();

        auto quadPoints = _quad->GetPoints();
        auto weights    = _quad->GetWeights();
        auto nq         = _quad->GetNq();

        // Iterate over boundaries
        for( std::map<int, Vector>::iterator it = _ghostCells.begin(); it != _ghostCells.end(); ++it ) {
            int idx_cell = it->first;    // Get Boundary cell index

            // Iterate over face cell faces

#pragma omp parallel for reduction( + : _curScalarOutflow )
            for( unsigned idx_nbr = 0; idx_nbr < neigbors[idx_cell].size(); ++idx_nbr ) {
                // Find face that points outward
                if( neigbors[idx_cell][idx_nbr] == nCells ) {
                    // Iterate over transport directions
                    for( unsigned idx_quad = 0; idx_quad < nq; ++idx_quad ) {
                        transportDirection =
                            normals[idx_cell][idx_nbr][0] * quadPoints[idx_quad][0] + normals[idx_cell][idx_nbr][1] * quadPoints[idx_quad][1];
                        // Find outward facing transport directions
                        if( transportDirection > 0.0 ) {
                            _curScalarOutflow += transportDirection * solution[idx_cell][idx_quad] * weights[idx_quad];    // Integrate flux
                        }
                    }
                }
            }
        }
    }
    // TODO define alternative for Moment solvers
}

void ProblemBase::ComputeTotalOutflow( double dT ) { _totalScalarOutflow += _curScalarOutflow * dT; }

void ProblemBase::ComputeMass( const Vector& scalarFlux ) {
    _mass = 0.0;

    auto areas      = _mesh->GetCellAreas();
    unsigned nCells = _mesh->GetNumCells();
#pragma omp parallel reduction( + : _mass )
    for( unsigned idx_cell = 0; idx_cell < nCells; ++idx_cell ) {
        _mass += scalarFlux[idx_cell] * areas[idx_cell];
    }
}

void ProblemBase::ComputeChangeRateFlux( const Vector& scalarFlux, const Vector& scalarFluxNew ) {
    _changeRateFlux = blaze::l2Norm( scalarFluxNew - scalarFlux );
}
