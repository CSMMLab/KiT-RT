#include "quadratures/qgausschebyshev1D.hpp"
#include "toolboxes/errormessages.hpp"

QGaussChebyshev1D::QGaussChebyshev1D( Config* settings ) : QuadratureBase( settings ) {
    SetName();
    CheckOrder();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
    _supportedDimensions = { 1 };
}

QGaussChebyshev1D::QGaussChebyshev1D( unsigned quadOrder ) : QuadratureBase( quadOrder ) {
    SetName();
    CheckOrder();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
}

void QGaussChebyshev1D::SetPointsAndWeights() {
    _pointsKarth.resize( _nq );
    _weights.resize( _nq );
    unsigned dim = 1;
    for( unsigned k = 1; k <= _nq; ++k ) {
        _pointsKarth[k - 1].resize( dim );
        _pointsKarth[k - 1][0] = std::cos( ( 2 * k - 1 ) * M_PI / ( 2 * _nq ) );
        _weights[k - 1]        = 1.0 / std::sqrt( 1 - _pointsKarth[k - 1][0] * _pointsKarth[k - 1][0] );
        std::cout << _pointsKarth[k - 1][0] << "\t" << _weights[k - 1] << std::endl;
    }
    _pointsSphere = _pointsKarth;
}

void QGaussChebyshev1D::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

bool QGaussChebyshev1D::CheckOrder() {
    return true;    // All orders viable
}

double QGaussChebyshev1D::Integrate( double ( *f )( double, double, double ) ) {    // Not Safe!
    double result = 0.0;
    double x      = 0.0;
    double y      = 0.0;
    double z      = 0.0;
    double w      = 0.0;
    for( unsigned i = 0; i < _nq; i++ ) {
        x = _pointsKarth[i][0];
        w = _weights[i];
        result += w * f( x, y, z );
    }
    return result;
}

double QGaussChebyshev1D::IntegrateSpherical( double ( *f )( double, double ) ) {
    ErrorMessages::Error( "This method is not applicable for 1D quadratures.\n", CURRENT_FUNCTION );
    return f( 0, 0 );
}

std::vector<double> QGaussChebyshev1D::Integrate( std::vector<double> ( *f )( double, double, double ), unsigned /* len */ ) {
    ErrorMessages::Error( "This method is not applicable for 1D quadratures.\n", CURRENT_FUNCTION );
    return { f( 0, 0, 0 ) };
}
