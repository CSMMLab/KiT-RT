#include "quadratures/qgausschebyshev1D.h"
#include "toolboxes/errormessages.h"

QGaussChebyshev1D::QGaussChebyshev1D( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    CheckOrder();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
}

QGaussChebyshev1D::QGaussChebyshev1D( unsigned quadOrder ) : QuadratureBase( quadOrder ) {
    SetName();
    CheckOrder();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
}

void QGaussChebyshev1D::SetPointsAndWeights() {
    _points.resize( _nq );
    _weights.resize( _nq );
    unsigned dim = 1;
    for( unsigned k = 1; k <= _nq; ++k ) {
        _points[k - 1].resize( dim );
        _points[k - 1][0] = std::cos( ( 2 * k - 1 ) * M_PI / ( 2 * _nq ) );
        _weights[k - 1]   = 1.0 / std::sqrt( 1 - _points[k - 1][0] * _points[k - 1][0] );
        std::cout << _points[k - 1][0] << "\t" << _weights[k - 1] << std::endl;
    }
    _pointsSphere = _points;
}

void QGaussChebyshev1D::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

bool QGaussChebyshev1D::CheckOrder() {
    return true;    // All orders viable
}

double QGaussChebyshev1D::Integrate( double( f )( double x0, double x1, double x2 ) ) {    // Not Safe!
    double result = 0.0;
    double x      = 0.0;
    double y      = 0.0;
    double z      = 0.0;
    double w      = 0.0;
    for( unsigned i = 0; i < _nq; i++ ) {
        x = _points[i][0];
        w = _weights[i];
        result += w * f( x, y, z );
    }
    return result;
}

double QGaussChebyshev1D::IntegrateSpherical( double( f )( double my, double phi ) ) {
    ErrorMessages::Error( "This method is not applicable for 1D quadratures.\n", CURRENT_FUNCTION );
    return f( 0, 0 );
}

std::vector<double> QGaussChebyshev1D::Integrate( std::vector<double>( f )( double x0, double x1, double x2 ), unsigned /* len */ ) {
    ErrorMessages::Error( "This method is not applicable for 1D quadratures.\n", CURRENT_FUNCTION );
    return { f( 0, 0, 0 ) };
}
