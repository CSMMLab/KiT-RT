#include "problems/isotropicsource2d.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"
#include "toolboxes/errormessages.hpp"

IsotropicSource2D::IsotropicSource2D( Config* settings, Mesh* mesh ) : ElectronRT( settings, mesh ) {}

IsotropicSource2D::~IsotropicSource2D() {}

std::vector<VectorVector> IsotropicSource2D::GetExternalSource( const Vector& energies ) {
    auto zeroVec = Vector( _settings->GetNQuadPoints(), 0.0 );
    auto uniform = std::vector<Vector>( _mesh->GetNumCells(), zeroVec );
    auto Q       = std::vector<VectorVector>( energies.size(), uniform );
    return Q;
}

VectorVector IsotropicSource2D::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids         = _mesh->GetCellMidPoints();
    double enterPositionX = 0.5;    // 0.0;
    double enterPositionY = 0.5;
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

std::vector<double> IsotropicSource2D::GetDensity( const VectorVector& /*cellMidPoints*/ ) {
    double rhoL = 1.0;
    double rhoR = 5.0;
    std::vector<double> rho( _settings->GetNCells(), rhoL );

    // use values rhoR on right third of domain
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        if( x >= 0.56 ) {
            rho[j] = rhoR;
        }
    }
    return rho;
}

IsotropicSource2D_Moment::IsotropicSource2D_Moment( Config* settings, Mesh* mesh ) : ElectronRT( settings, mesh ) {}

IsotropicSource2D_Moment::~IsotropicSource2D_Moment() {}

std::vector<VectorVector> IsotropicSource2D_Moment::GetExternalSource( const Vector& energies ) {
    auto zeroVec = Vector( _settings->GetNQuadPoints(), 0.0 );
    auto uniform = std::vector<Vector>( _mesh->GetNumCells(), zeroVec );
    auto Q       = std::vector<VectorVector>( energies.size(), uniform );
    // ErrorMessages::Error( "Function not yet implemented.", CURRENT_FUNCTION );
    return Q;
}

VectorVector IsotropicSource2D_Moment::SetupIC() {
    VectorVector psi( _mesh->GetNumCells(), Vector( _settings->GetNQuadPoints(), 1e-10 ) );
    auto cellMids         = _mesh->GetCellMidPoints();
    double enterPositionX = 0.5;    // 0.0;
    double enterPositionY = 0.5;

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

std::vector<double> IsotropicSource2D_Moment::GetDensity( const VectorVector& /*cellMidPoints*/ ) {
    double rhoL = 1.0;
    double rhoR = 5.0;
    std::vector<double> rho( _settings->GetNCells(), rhoL );

    // use values rhoR on right third of domain
    auto cellMids = _mesh->GetCellMidPoints();
    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        if( x >= 0.56 ) {
            rho[j] = rhoR;
        }
    }
    return rho;
}
