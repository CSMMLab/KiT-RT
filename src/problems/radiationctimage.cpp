#include "problems/radiationctimage.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include "common/mesh.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/interpolation.hpp"
#include "velocitybasis/sphericalbase.hpp"
#include "quadratures/quadraturebase.hpp"
#include "velocitybasis/sphericalbase.hpp"
#include "velocitybasis/sphericalharmonics.hpp"
#include "quadratures/qgausslegendretensorized.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"
#include "solvers/csdpn_starmap_constants.hpp"

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
    VectorVector cellMids         = _mesh->GetCellMidPoints();
    double s                      = 0.1;
    double enterPositionX = 2.5;    // 0.0;
    double enterPositionY = 5.8;
    double meanDir = M_PI/2;
    double s_ang = 0.1;
    double epsilon = 1e-3;
    QuadratureBase* quad          = QuadratureBase::Create( _settings );
    VectorVector quadPointsSphere = quad->GetPointsSphere();

    for( unsigned j = 0; j < cellMids.size(); ++j ) {
            double x = cellMids[j][0] - enterPositionX;
            double y = cellMids[j][1] - enterPositionY;
        // anisotropic inflow that concentrates all particles on the last quadrature point
    //    for( unsigned idx_quad = 0; idx_quad < _settings->GetNQuadPoints(); idx_quad++ ) {
    //           if( quadPointsSphere[idx_quad][1] > 4*M_PI/3 && quadPointsSphere[idx_quad][1] < 5*M_PI/3) {    // if my >0 
    //                 psi[j][idx_quad] = std::max( 1.0 / ( s * sqrt( 2 * M_PI ) )  * std::exp( -( x * x + y * y ) / ( 2 * s * s ) ), epsilon );
    //         }
    //         }
 // normal distribution also in angle
            
        for( unsigned idx_quad = 0; idx_quad < _settings->GetNQuadPoints(); idx_quad++ ) {
                double ang = quadPointsSphere[idx_quad][1] - meanDir;
                psi[j][idx_quad]= std::max( 1.0 / ( s * s_ang * sqrt( 8 * M_PI * M_PI * M_PI))  * std::exp( -( x * x + y * y ) / ( 2 * s * s ) )* std::exp( -( ang * ang ) / (2* s_ang ) ), epsilon );
        }
    }


    delete quad;
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
        result[i] = std::clamp( interp( cellMidPoints[i][0], cellMidPoints[i][1] )*1.85, 0.05, 1.85 ); //Scale densities for CT to be between 0 (air) and 1.85 (bone)
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
      if( _settings->GetSphericalBasisName() == SPHERICAL_HARMONICS  ) {
        // In case of PN, spherical basis is per default SPHERICAL_HARMONICS in 3 velocity dimensions
        SphericalBase* tempBase  = new SphericalHarmonics( _settings->GetMaxMomentDegree(), 3 );
        unsigned ntotalEquations = tempBase->GetBasisSize();

        double epsilon = 1e-3;

        VectorVector initialSolution( _mesh->GetNumCells(), Vector( ntotalEquations, 0) );    
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
        double enterPositionX = 2.5;    // 0.0;
        double enterPositionY = 5.8;
        double meanDir = M_PI/2;
        double s_ang = 0.1;

        for( unsigned idx_cell = 0; idx_cell < cellMids.size(); ++idx_cell ) {
            double x = cellMids[idx_cell][0] - enterPositionX;
            double y = cellMids[idx_cell][1] - enterPositionY;
            // anisotropic, forward-directed particle inflow
            // for( unsigned idx_quad = 0; idx_quad < _settings->GetNQuadPoints(); idx_quad++ ) {
            //     if( quadPointsSphere[idx_quad][1] > M_PI/3 && quadPointsSphere[idx_quad][1] < 2*M_PI/3) {    // if my >0 
            //         cellKineticDensity[idx_quad] = std::max( 1.0 / ( s * sqrt( 2 * M_PI ) )  * std::exp( -( x * x + y * y ) / ( 2 * s * s ) ), epsilon );
            //     }
            //  }

            // normal distribution also in angle
           
            for( unsigned idx_quad = 0; idx_quad < _settings->GetNQuadPoints(); idx_quad++ ) {
                    double ang = quadPointsSphere[idx_quad][1] - meanDir;
                    cellKineticDensity[idx_quad] = std::max( 1.0 / ( s * s *s_ang * sqrt( 8 * M_PI * M_PI * M_PI))  * std::exp( -( x * x + y * y ) / ( 2 * s * s ) )* std::exp( -( ang * ang ) / (2* s_ang ) ), epsilon );
            }
            
            // Compute moments of this kinetic density
            for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
                initialSolution[idx_cell] += cellKineticDensity[idx_quad] * w[idx_quad] * moments[idx_quad];
            }
        }
        //TextProcessingToolbox::PrintVectorVector(initialSolution);
        //exit(0);
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
        double enterPositionX = 2.5;    // 0.0;
        double enterPositionY = 5.8;
        double meanDir = M_PI/2;
        double s_ang = M_PI *3;


        for( unsigned idx_cell = 0; idx_cell < cellMids.size(); ++idx_cell ) {
            double x = cellMids[idx_cell][0] - enterPositionX;
            double y = cellMids[idx_cell][1] - enterPositionY;
            
            // anisotropic, forward-directed particle inflow
            for( unsigned idx_quad = 0; idx_quad < _settings->GetNQuadPoints(); idx_quad++ ) {
                if( quadPointsSphere[idx_quad][1] > M_PI/3 && quadPointsSphere[idx_quad][1] < 2*M_PI/3) {    // if my >0 
                    cellKineticDensity[idx_quad] = std::max( 1.0 / ( s * sqrt( 2 * M_PI ) )  * std::exp( -( x * x + y * y ) / ( 2 * s * s ) ), epsilon );
                }
            }
            // normal distribution also in angle
           
            // for( unsigned idx_quad = 0; idx_quad < _settings->GetNQuadPoints(); idx_quad++ ) {
            //         double ang = quadPointsSphere[idx_quad][1] - meanDir;
            //         cellKineticDensity[idx_quad] = std::max( 1.0 / ( s * s*  s_ang * sqrt(  8 * M_PI * M_PI * M_PI ) )  * std::exp( -( x * x + y * y ) / ( 2 * s * s ) )* std::exp( -( ang * ang ) / (2* s_ang ) ), epsilon );
            // }
            // Compute moments of this kinetic density
            for( unsigned idx_quad = 0; idx_quad < quad->GetNq(); idx_quad++ ) {
                initialSolution[idx_cell] += cellKineticDensity[idx_quad] * w[idx_quad] * moments[idx_quad];
            }
        }
        delete quad;
        return initialSolution;
    }
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
        result[i] = std::clamp( interp( cellMidPoints[i][0], cellMidPoints[i][1] )*1.85, 0.05, 1.85 ); //Scale densities for CT to be between 0 (air) and 1.85 (bone)
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
double RadiationCTImage::NormPDF( double x, double mu, double sigma ) {
    return INV_SQRT_2PI / sigma * std::exp( -( ( x - mu ) * ( x - mu ) ) / ( 2.0 * sigma * sigma ) );
}