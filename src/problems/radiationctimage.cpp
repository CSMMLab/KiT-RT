#include "problems/radiationctimage.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/interpolation.hpp"
#include "velocitybasis/sphericalbase.hpp"

#include <fstream>

RadiationCTImage::RadiationCTImage( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {}

RadiationCTImage::~RadiationCTImage() {}

std::vector<VectorVector> RadiationCTImage::GetExternalSource( const Vector& energies ) {
    auto zeroVec = Vector( _settings->GetNQuadPoints(), 0.0 );
    auto uniform = std::vector<Vector>( _mesh->GetNumCells(), zeroVec );
    auto Q       = std::vector<VectorVector>( energies.size(), uniform );
    return Q;
}

VectorVector RadiationCTImage::SetupIC() {
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

std::vector<double> RadiationCTImage::GetDensity( const VectorVector& /*cellMidPoints*/ ) {
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
        result[i] = std::clamp( interp( cellMidPoints[i][0], cellMidPoints[i][1] )*1.85, 0.4, 1.85 ); //Scale densities for CT to be between 0 (air) and 1.85 (bone)
    }
    return result;
}
VectorVector RadiationCTImage::GetScatteringXS( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( 1, Vector( 1, 0.0 ) );
}

VectorVector RadiationCTImage::GetTotalXS( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( 1, Vector( 1, 0.0 ) );
}
RadiationCTImage_Moment::RadiationCTImage_Moment( Config* settings, Mesh* mesh ) : ProblemBase( settings, mesh ) {}

RadiationCTImage_Moment::~RadiationCTImage_Moment() {}

std::vector<VectorVector> RadiationCTImage_Moment::GetExternalSource( const Vector& energies ) {
    SphericalBase* tempBase  = SphericalBase::Create( _settings );
    unsigned ntotalEquations = tempBase->GetBasisSize();
    delete tempBase;

    Vector zeroVec              = Vector( ntotalEquations, 0.0 );
    VectorVector uniform        = VectorVector( _mesh->GetNumCells(), zeroVec );
    std::vector<VectorVector> Q = std::vector<VectorVector>( energies.size(), uniform );
    return Q;
}

VectorVector RadiationCTImage_Moment::SetupIC() {
    // Compute number of equations in the system

    // In case of PN, spherical basis is per default SPHERICAL_HARMONICS
    SphericalBase* tempBase  = SphericalBase::Create( _settings );
    unsigned ntotalEquations = tempBase->GetBasisSize();

    VectorVector initial_sol( _mesh->GetNumCells(), Vector( ntotalEquations, 0 ) );    // zero could lead to problems?
    VectorVector cellMids = _mesh->GetCellMidPoints();

    Vector uIC( ntotalEquations, 0 );

    if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
        QuadratureBase* quad          = QuadratureBase::Create( _settings );
        VectorVector quadPointsSphere = quad->GetPointsSphere();
        Vector w                      = quad->GetWeights();

        double my, phi;
        VectorVector moments = VectorVector( quad->GetNq(), Vector( tempBase->GetBasisSize(), 0.0 ) );

        for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
            my                = quadPointsSphere[idx_quad][0];
            phi               = quadPointsSphere[idx_quad][1];
            moments[idx_quad] = tempBase->ComputeSphericalBasis( my, phi );
        }
        // Integrate <1*m> to get factors for monomial basis in isotropic scattering
        for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
            uIC += w[idx_quad] * moments[idx_quad];
        }
        delete quad;
    }

    // Initial condition is dirac impulse at (x,y) = (0,0) ==> constant in angle ==> all moments - exept first - are zero.
    double t       = 3.2e-4;    // pseudo time for gaussian smoothing (Approx to dirac impulse)
    double epsilon = 1e-4;      // minimal value for first moment to avoid div by zero error

    for( unsigned j = 0; j < cellMids.size(); ++j ) {
        double x = cellMids[j][0];
        double y = cellMids[j][1];    // (x- 0.5) * (x- 0.5)

        double kinetic_density = std::max( 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( x * x + y * y ) / ( 4 * t ) ), epsilon );

        if( _settings->GetSphericalBasisName() == SPHERICAL_MONOMIALS ) {
            initial_sol[j] = kinetic_density * uIC / uIC[0];    // Remember scaling
        }
        if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS ) {
            initial_sol[j][0] = kinetic_density;
        }
    }
    delete tempBase;    // Only temporally needed
    std::cout << "not correct iC. but gets overwritten in csdpn_jl constructor\n";
    return initial_sol;
}

std::vector<double> RadiationCTImage_Moment::GetDensity( const VectorVector& /*cellMidPoints*/ ) {
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
        result[i] = std::clamp( interp( cellMidPoints[i][0], cellMidPoints[i][1] )*1.85, 0.4, 1.85 ); //Scale densities for CT to be between 0 (air) and 1.85 (bone)
    }
    return result;
}
VectorVector RadiationCTImage_Moment::GetScatteringXS( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( 1, Vector( 1, 0.0 ) );
}

VectorVector RadiationCTImage_Moment::GetTotalXS( const Vector& /*energies*/ ) {
    // @TODO
    // Specified in subclasses
    return VectorVector( 1, Vector( 1, 0.0 ) );
}
