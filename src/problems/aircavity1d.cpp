#include "problems/aircavity1d.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"
#include "quadratures/qgausslegendretensorized.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/sphericalbase.hpp"
#include "toolboxes/sphericalharmonics.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

AirCavity1D::AirCavity1D( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _sigmaS = settings->GetSigmaS(); }

AirCavity1D::~AirCavity1D() {}

std::vector<VectorVector> AirCavity1D::GetExternalSource( const Vector& energies ) {
    return std::vector<VectorVector>( energies.size(), std::vector<Vector>( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 0.0 ) ) );
}

VectorVector AirCavity1D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    VectorVector cellMids         = _mesh->GetCellMidPoints();
    double s                      = 0.1;
    QuadratureBase* quad          = QuadratureBase::Create( _settings );
    VectorVector quadPointsSphere = quad->GetPointsSphere();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        // anisotropic inflow that concentrates all particles on the last quadrature point
        for( unsigned idx_quad = 0; idx_quad < _settings->GetNQuadPoints(); idx_quad++ ) {
            if( quadPointsSphere[idx_quad][0] > 0.5 ) {    // if my >0
                psi[j][idx_quad] = 1.0 / ( s * sqrt( 2 * M_PI ) ) * std::exp( -x * x / ( 2 * s * s ) );
            }
        }
    }
    delete quad;
    return psi;
}

std::vector<double> AirCavity1D::GetDensity( const VectorVector& cellMidPoints ) {
    std::vector<double> densities( cellMidPoints.size(), 1.0 );
    for( unsigned j = 0; j < cellMidPoints.size(); ++j ) {
        if( cellMidPoints[j][0] > 1.5 - 2.5 && cellMidPoints[j][0] < 2.0 - 2.5 ) densities[j] = 0.01;
    }
    return densities;
}

VectorVector AirCavity1D::GetScatteringXS( const Vector& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) );
}

VectorVector AirCavity1D::GetTotalXS( const Vector& energies ) { return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) ); }

// ------ Moment version ---

AirCavity1D_Moment::AirCavity1D_Moment( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) { _sigmaS = settings->GetSigmaS(); }

AirCavity1D_Moment::~AirCavity1D_Moment() {}

std::vector<VectorVector> AirCavity1D_Moment::GetExternalSource( const Vector& /*energies*/ ) {
    SphericalBase* tempBase  = SphericalBase::Create( _settings );
    unsigned ntotalEquations = tempBase->GetBasisSize();
    VectorVector Q( _mesh->GetNumCells(), Vector( ntotalEquations, 0.0 ) );    // zero could lead to problems?
    delete tempBase;                                                           // Only temporally needed
    return std::vector<VectorVector>( 1u, Q );
}

VectorVector AirCavity1D_Moment::SetupIC() {
    if( _settings->GetSolverName() == PN_SOLVER ) {
        // In case of PN, spherical basis is per default SPHERICAL_HARMONICS in 3 velocity dimensions

        SphericalBase* tempBase  = new SphericalHarmonics( _settings->GetMaxMomentDegree(), 3 );
        unsigned ntotalEquations = tempBase->GetBasisSize();

        double epsilon = 1e-3;

        VectorVector initialSolution( _mesh->GetNumCells(), Vector( ntotalEquations, 0.0 ) );    // zero could lead to problems?
        VectorVector cellMids = _mesh->GetCellMidPoints();

        QuadratureBase* quad          = new QGaussLegendreTensorized( _settings );
        VectorVector quadPointsSphere = quad->GetPointsSphere();
        Vector w                      = quad->GetWeights();
        Vector cellKineticDensity( quad->GetNq(), epsilon );

        // compute moment basis
        VectorVector moments = VectorVector( quad->GetNq(), Vector( ntotalEquations, 0.0 ) );
        double my, phi;    // quadpoints in spherical coordinates
        for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
            my                = quadPointsSphere[idx_quad][0];
            phi               = quadPointsSphere[idx_quad][1];
            moments[idx_quad] = tempBase->ComputeSphericalBasis( my, phi );
        }
        delete tempBase;
        double s = 0.1;
        // create kinetic density( SN initial condition )
        for( unsigned idx_cell = 0; idx_cell < cellMids.size(); ++idx_cell ) {
            double x = cellMids[idx_cell][0];
            // anisotropic inflow that concentrates all particles on the last quadrature point
            for( unsigned idx_quad = 0; idx_quad < _settings->GetNQuadPoints(); idx_quad++ ) {
                if( quadPointsSphere[idx_quad][0] > 0.5 ) {    // if my >0
                    cellKineticDensity[idx_quad] = std::max( 1.0 / ( s * sqrt( 2 * M_PI ) ) * std::exp( -x * x / ( 2 * s * s ) ), epsilon );
                }
            }
            // Compute moments of this kinetic density
            for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
                initialSolution[idx_cell] += cellKineticDensity[idx_quad] * w[idx_quad] * moments[idx_quad];
            }
        }
        delete quad;
        return initialSolution;
    }
    else {
        SphericalBase* tempBase  = SphericalBase::Create( _settings );
        unsigned ntotalEquations = tempBase->GetBasisSize();

        double epsilon = 1e-3;

        Vector cellKineticDensity( _settings->GetNQuadPoints(), epsilon );
        VectorVector initialSolution( _mesh->GetNumCells(), Vector( ntotalEquations, 0 ) );    // zero could lead to problems?
        VectorVector cellMids = _mesh->GetCellMidPoints();

        QuadratureBase* quad          = QuadratureBase::Create( _settings );
        VectorVector quadPointsSphere = quad->GetPointsSphere();
        Vector w                      = quad->GetWeights();

        // compute moment basis
        VectorVector moments = VectorVector( quad->GetNq(), Vector( ntotalEquations, 0.0 ) );
        double my, phi;    // quadpoints in spherical coordinates
        for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
            my                = quadPointsSphere[idx_quad][0];
            phi               = quadPointsSphere[idx_quad][1];
            moments[idx_quad] = tempBase->ComputeSphericalBasis( my, phi );
        }
        delete tempBase;

        double s = 0.1;
        // create kinetic density( SN initial condition )
        for( unsigned idx_cell = 0; idx_cell < cellMids.size(); ++idx_cell ) {
            double x = cellMids[idx_cell][0];
            // anisotropic inflow that concentrates all particles on the last quadrature point
            for( unsigned idx_quad = 0; idx_quad < _settings->GetNQuadPoints(); idx_quad++ ) {
                if( quadPointsSphere[idx_quad][0] > 0.5 ) {    // if my >0
                    cellKineticDensity[idx_quad] = std::max( 1.0 / ( s * sqrt( 2 * M_PI ) ) * std::exp( -x * x / ( 2 * s * s ) ), epsilon );
                }
            }
            // Compute moments of this kinetic density
            for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
                initialSolution[idx_cell] += cellKineticDensity[idx_quad] * w[idx_quad] * moments[idx_quad];
            }
        }
        delete quad;
        return initialSolution;
    }
}

std::vector<double> AirCavity1D_Moment::GetDensity( const VectorVector& cellMidPoints ) {
    std::vector<double> densities( cellMidPoints.size(), 1.0 );

    for( unsigned j = 0; j < cellMidPoints.size(); ++j ) {
        if( cellMidPoints[j][0] > 1.5 - 2.5 && cellMidPoints[j][0] < 2.0 - 2.5 ) densities[j] = 0.01;
    }
    return densities;
}

VectorVector AirCavity1D_Moment::GetScatteringXS( const Vector& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) );
}

VectorVector AirCavity1D_Moment::GetTotalXS( const Vector& energies ) {
    return VectorVector( energies.size(), Vector( _mesh->GetNumCells(), _sigmaS ) );
}
