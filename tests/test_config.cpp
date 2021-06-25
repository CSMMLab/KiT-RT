#include "catch.hpp"
#include "common/config.h"

TEST_CASE( "Read in Config Template" ) {
    std::string filename = std::string( TESTS_PATH ) + "input/unit_tests/common/unit_config.cfg";

    // Load Settings from File
    Config* config = new Config( filename );

    // Test all set configurations
    bool allSatisfied = true;

    if( config->GetCFL() != 0.4 ) allSatisfied = false;
    if( config->GetTEnd() != 0.3 ) allSatisfied = false;

    if( config->GetProblemName() != PROBLEM_LineSource ) allSatisfied = false;
    if( config->GetKernelName() != KERNEL_Isotropic ) allSatisfied = false;
    if( config->GetQuadName() != QUAD_MonteCarlo ) allSatisfied = false;

    if( config->GetQuadOrder() != 5000 ) allSatisfied = false;

    // Check boundary conditions
    if( config->GetBoundaryType( "DirichletTestMarker2" ) != DIRICHLET ) allSatisfied = false;
    if( config->GetBoundaryType( "DirichletTestMarker1" ) != DIRICHLET ) allSatisfied = false;
    if( config->GetBoundaryType( "NeumannTestMarker1" ) != NEUMANN ) allSatisfied = false;
    if( config->GetBoundaryType( "NeumannTestMarker2" ) != NEUMANN ) allSatisfied = false;

    REQUIRE( allSatisfied );
}
