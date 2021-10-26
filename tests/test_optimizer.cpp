#include <numeric>

#include "catch.hpp"
#include "common/config.hpp"
#include "optimizers/optimizerbase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "toolboxes/sphericalharmonics.hpp"

TEST_CASE( "Test the Newton Optimizer", "[optimizers]" ) {
    std::string filename = std::string( TESTS_PATH ) + "input/unit_tests/optimizers/unit_optimizerNewton.cfg";

    // Load Settings from File
    Config* config = new Config( filename );

    // Get Basis
    SphericalHarmonics basis( config->GetMaxMomentDegree() );

    // Get Quadrature
    QuadratureBase* quad = QuadratureBase::Create( config );

    // Get Optimizer (Newton)
    OptimizerBase* optimizer = OptimizerBase::Create( config );

    // Get dummy Moment Vector
    unsigned nTotalEntries = basis.GetGlobalIndexBasis( config->GetMaxMomentDegree(), config->GetMaxMomentDegree() ) + 1;    // = 4
    Vector u( nTotalEntries, -1.5 );
    u[1] = 0.0;
    u[2] = 1.0;

    // Get inital guess for solution
    Vector alpha( nTotalEntries, 27.0 );

    // Get Moments

    VectorVector moments = VectorVector( quad->GetNq() );
    double my, phi;
    VectorVector quadPointsSphere = quad->GetPointsSphere();
    for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
        my  = quadPointsSphere[idx_quad][0];
        phi = quadPointsSphere[idx_quad][1];

        moments[idx_quad] = basis.ComputeSphericalBasis( my, phi );
    }

    // Solve
    optimizer->Solve( alpha, u, moments );

    REQUIRE( std::fabs( norm( alpha - u ) ) < config->GetNewtonOptimizerEpsilon() );    // alpha = u for quadratic entropy
}
