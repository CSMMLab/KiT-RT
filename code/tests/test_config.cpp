#include "catch.hpp"
#include "settings/config.h"
#include "settings/globalconstants.h"

TEST_CASE( "Read in Config Template" ) {

    std::string filename = "input/configTest.cfg";
    char config_file_name[MAX_STRING_SIZE];

    /*--- Load in the number of zones and spatial dimensions in the mesh file (If no config
     file is specified, default.cfg is used) ---*/
    strcpy( config_file_name, filename.c_str() );

    // Load Settings from File
    Config* config = new Config( config_file_name );

    // Test all set configurations
    bool allSatisfied = true;

    if( config->GetMeshFile().compare( "testMesh.su2" ) != 0 ) allSatisfied = false;
    if( config->GetOutputDir().compare( "../result" ) != 0 ) allSatisfied = false;
    if( config->GetOutputFile().compare( "testOutput" ) != 0 ) allSatisfied = false;
    if( config->GetLogDir().compare( "../result/logs" ) != 0 ) allSatisfied = false;

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
