#include "problems/isotropicsource2d_ct.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/interpolation.hpp"
#include <fstream>
#include <numeric>

IsotropicSource2D_CT::IsotropicSource2D_CT( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {}

IsotropicSource2D_CT::~IsotropicSource2D_CT() {}

std::vector<VectorVector> IsotropicSource2D_CT::GetExternalSource( const Vector& energies ) {
    auto zeroVec = Vector( _settings->GetNQuadPoints(), 0.0 );
    auto uniform = std::vector<Vector>( _mesh->GetNumCells(), zeroVec );
    auto Q       = std::vector<VectorVector>( energies.size(), uniform );
    return Q;
}

VectorVector IsotropicSource2D_CT::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids         = _mesh->GetCellMidPoints();
    double enterPositionX = 0.5 * 10;    // 0.0;
    double enterPositionY = 0.5 * 10;
    // auto boundaryCells    = _mesh->GetBoundaryTypes();
    // Case 1: Ingoing radiation in just one cell
    // find cell that best matches enter position
    // double dist = 1000.0;
    // unsigned indexSource = 0;
    /*
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        // if( boundaryCells[j] == BOUNDARY_TYPE::DIRICHLET ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];
        if( x >= 0.49 && x <= 0.5 && y >= 0.49 && y <= 0.5 ) {
            psi[j] = Vector( _settings->GetNQuadPoints(), 1.0 );
        }
        //}
    }*/
    // psi[indexSource] = Vector( _settings->GetNQuadPoints(), 1.0 );

    // Case 2: Ingoing radiation as Gauss curve
    double t = 1e-5;    // pseudo time for gaussian smoothing
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0] - enterPositionX;
        double y = cellMids[j][1] - enterPositionY;
        psi[j]   = Vector( _settings->GetNQuadPoints(), 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) ) );
    }
    return psi;
}

std::vector<double> IsotropicSource2D_CT::GetDensity( const VectorVector& /*cellMidPoints*/ ) {
    std::string imageFile = _settings->GetCTFile();
    std::string meshFile  = _settings->GetMeshFile();
    Matrix gsImage        = createSU2MeshFromImage( imageFile, meshFile );
    auto bounds           = _mesh->GetBounds();
    auto cellMidPoints    = _mesh->GetCellMidPoints();

    double xMin = bounds[0].first;
    double xMax = bounds[0].second;
    double yMin = bounds[1].first;
    double yMax = bounds[1].second;

    unsigned m = gsImage.rows();
    std::cout << "Number rows: " << m << std::endl;
    unsigned n = gsImage.columns();
    std::cout << "Number columns: " << n << std::endl;

    Vector x( m ), y( n );
    for( unsigned i = 0; i < m; ++i ) {
        x[i] = xMin + static_cast<double>( i ) / static_cast<double>( m - 1 ) * ( xMax - xMin );
    }
    for( unsigned i = 0; i < n; ++i ) y[i] = yMin + static_cast<double>( i ) / static_cast<double>( n - 1 ) * ( yMax - yMin );

    Interpolation interp( x, y, gsImage );
    std::vector<double> result( _mesh->GetNumCells(), 0.0 );
    std::ofstream fout;
    fout.open( "density_test.txt" );
    std::ofstream fout1;
    fout1.open( "x_test.txt" );
    std::ofstream fout2;
    fout2.open( "y_test.txt" );
    for( unsigned i = 0; i < _mesh->GetNumCells(); ++i ) {
        result[i] = std::clamp( interp( cellMidPoints[i][0], cellMidPoints[i][1] ), 0.6, 1.85 );
        fout1 << cellMidPoints[i][0] << std::endl;
        fout2 << cellMidPoints[i][1] << std::endl;
        fout << result[i] << std::endl;
    }
    fout.close();
    return result;
}

IsotropicSource2D_CT_Moment::IsotropicSource2D_CT_Moment( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {}

IsotropicSource2D_CT_Moment::~IsotropicSource2D_CT_Moment() {}

std::vector<VectorVector> IsotropicSource2D_CT_Moment::GetExternalSource( const Vector& energies ) {
    auto zeroVec = Vector( _settings->GetNQuadPoints(), 0.0 );
    auto uniform = std::vector<Vector>( _mesh->GetNumCells(), zeroVec );
    auto Q       = std::vector<VectorVector>( energies.size(), uniform );
    // ErrorMessages::Error( "Function not yet implemented.", CURRENT_FUNCTION );
    return Q;
}

VectorVector IsotropicSource2D_CT_Moment::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids         = _mesh->GetCellMidPoints();
    double enterPositionX = 0.5 * 10;    // 0.0;
    double enterPositionY = 0.5 * 10;

    // Case 2: Ingoing radiation as Gauss curve
    double t = 1e-5;    // pseudo time for gaussian smoothing
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0] - enterPositionX;
        double y = cellMids[j][1] - enterPositionY;
        psi[j]   = Vector( _settings->GetNQuadPoints(), 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) ) );
    }
    std::cout << "not correct iC. but gets overwritten in csdpn_jl constructor\n";
    return psi;
}

std::vector<double> IsotropicSource2D_CT_Moment::GetDensity( const VectorVector& /*cellMidPoints*/ ) {
    std::string imageFile = _settings->GetCTFile();
    std::string meshFile  = _settings->GetMeshFile();
    Matrix gsImage        = createSU2MeshFromImage( imageFile, meshFile );
    auto bounds           = _mesh->GetBounds();
    auto cellMidPoints    = _mesh->GetCellMidPoints();

    double xMin = bounds[0].first;
    double xMax = bounds[0].second;
    double yMin = bounds[1].first;
    double yMax = bounds[1].second;

    unsigned m = gsImage.rows();
    unsigned n = gsImage.columns();

    Vector x( m ), y( n );
    for( unsigned i = 0; i < m; ++i ) {
        x[i] = xMin + static_cast<double>( i ) / static_cast<double>( m - 1 ) * ( xMax - xMin );
    }
    for( unsigned i = 0; i < n; ++i ) y[i] = yMin + static_cast<double>( i ) / static_cast<double>( n - 1 ) * ( yMax - yMin );

    Interpolation interp( x, y, gsImage );
    std::vector<double> result( _mesh->GetNumCells(), 0.0 );
    for( unsigned i = 0; i < _mesh->GetNumCells(); ++i ) {
        result[i] = std::clamp( interp( cellMidPoints[i][0], cellMidPoints[i][1] ), 0.4, 1.85 );
    }
    std::cout << "**** Test 2 ****";
    return result;
}
