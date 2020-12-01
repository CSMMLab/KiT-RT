#include "quadratures/quadraturebase.h"
#include "common/config.h"
#include "quadratures/qgausslegendre1D.h"
#include "quadratures/qgausslegendretensorized.h"
#include "quadratures/qldfesa.h"
#include "quadratures/qlebedev.h"
#include "quadratures/qlevelsymmetric.h"
#include "quadratures/qmontecarlo.h"
#include "quadratures/qproduct.h"
#include "toolboxes/errormessages.h"

QuadratureBase::QuadratureBase( Config* settings ) {
    _settings = settings;
    _order    = settings->GetQuadOrder();
}

QuadratureBase::QuadratureBase( unsigned order ) : _order( order ) { _settings = nullptr; }

QuadratureBase* QuadratureBase::CreateQuadrature( Config* settings ) {
    QUAD_NAME name = settings->GetQuadName();

    switch( name ) {
        case QUAD_MonteCarlo: return new QMonteCarlo( settings );
        case QUAD_GaussLegendreTensorized: return new QGaussLegendreTensorized( settings );
        case QUAD_GaussLegendre1D: return new QGaussLegendre1D( settings );
        case QUAD_LevelSymmetric: return new QLevelSymmetric( settings );
        case QUAD_LDFESA: return new QLDFESA( settings );
        case QUAD_Lebedev: return new QLebedev( settings );
        case QUAD_Product: return new ProductQuadrature( settings );
        default: ErrorMessages::Error( "Creator for the chose quadrature does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    }
}

QuadratureBase* QuadratureBase::CreateQuadrature( QUAD_NAME name, unsigned quadOrder ) {

    switch( name ) {
        case QUAD_MonteCarlo: return new QMonteCarlo( quadOrder );
        case QUAD_GaussLegendreTensorized:
            ErrorMessages::Error( "This quadrature must be initialized with a settings constructor!", CURRENT_FUNCTION );
            break;
        case QUAD_GaussLegendre1D: return new QGaussLegendre1D( quadOrder );
        case QUAD_LevelSymmetric: return new QLevelSymmetric( quadOrder );
        case QUAD_LDFESA: return new QLDFESA( quadOrder );
        case QUAD_Lebedev: return new QLebedev( quadOrder );
        case QUAD_Product: return new ProductQuadrature( quadOrder );
        default: ErrorMessages::Error( "Creator for the chose quadrature does not yet exist. This is is the fault of the coder!", CURRENT_FUNCTION );
    }
}

double QuadratureBase::Integrate( double( f )( double x0, double x1, double x2 ) ) {
    double result = 0.0;
    for( unsigned i = 0; i < _nq; i++ ) {
        double x = _points[i][0];
        double y = _points[i][1];
        double z = _points[i][2];
        double w = _weights[i];
        result += w * f( x, y, z );
    }
    return result;
}

VectorVector QuadratureBase::GetPointsSphere() const {
    ErrorMessages::Error( "Quadrature points in spherical coordinates are not supported\nfor this quadrature. Exiting", CURRENT_FUNCTION );
    return _pointsSphere;
}

std::vector<double> QuadratureBase::Integrate( std::vector<double>( f )( double x0, double x1, double x2 ), unsigned len ) {
    std::vector<double> result( len, 0.0 );
    std::vector<double> funcEval( len, 0.0 );

    for( unsigned i = 0; i < _nq; i++ ) {
        double x = _points[i][0];
        double y = _points[i][1];
        double z = _points[i][2];
        double w = _weights[i];
        funcEval = f( x, y, z );
        for( unsigned idx_len = 0; idx_len < len; idx_len++ ) {
            result[idx_len] += w * funcEval[idx_len];
        }
    }
    return result;
}

double QuadratureBase::SumUpWeights() { return sum( _weights ); }

void QuadratureBase::PrintWeights() {
    auto log = spdlog::get( "event" );
    for( unsigned i = 0; i < _nq; i++ ) {
        double w = _weights[i];
        log->info( w );
    }
}

void QuadratureBase::PrintPoints() {
    auto log = spdlog::get( "event" );
    for( unsigned i = 0; i < _nq; i++ ) {
        double x = _points[i][0];
        double y = _points[i][1];
        double z = _points[i][2];
        log->info( "{0}, {1}, {2}", x, y, z );
    }
}
void QuadratureBase::PrintPointsAndWeights() {
    auto log = spdlog::get( "event" );
    for( unsigned i = 0; i < _nq; i++ ) {
        double x = _points[i][0];
        double y = _points[i][1];
        double z = _points[i][2];
        double w = _weights[i];
        log->info( "{0}, {1}, {2}, {3}", x, y, z, w );
    }
}
