#include "common/config.hpp"
#include "quadratures/qgausslegendretensorized.hpp"
#include "toolboxes/errormessages.hpp"

QGaussLegendreTensorized2D::QGaussLegendreTensorized2D( Config* settings ) : QGaussLegendreTensorized( settings ) {
    SetName();
    CheckOrder();
    SetNq();
    SetPointsAndWeights();
    _supportedDimensions = { 2 };
}

void QGaussLegendreTensorized2D::SetNq() { _nq = pow( GetOrder(), 2 ); }

bool QGaussLegendreTensorized2D::CheckOrder() {
    if( _order % 2 == 1 ) {    // order needs to be even
        ErrorMessages::Error( "ERROR! Order " + std::to_string( _order ) + " for " + GetName() + " not available. \n Order must be an even number. ",
                              CURRENT_FUNCTION );
    }
    return true;
}

void QGaussLegendreTensorized2D::SetPointsAndWeights() {
    Vector nodes1D( _order ), weights1D( _order );

    // construct companion matrix
    Matrix CM( _order, _order, 0.0 );
    for( unsigned i = 0; i < _order - 1; ++i ) {
        CM( i + 1, i ) = std::sqrt( 1 / ( 4 - 1 / std::pow( static_cast<double>( i + 1 ), 2 ) ) );
        CM( i, i + 1 ) = std::sqrt( 1 / ( 4 - 1 / std::pow( static_cast<double>( i + 1 ), 2 ) ) );
    }

    // compute eigenvalues and -vectors of the companion matrix
    auto evSys = ComputeEigenValTriDiagMatrix( CM );

    for( unsigned i = 0; i < _order; ++i ) {
        if( std::fabs( evSys.first[i] ) < 1e-15 )    // avoid rounding errors
            nodes1D[i] = 0;
        else
            nodes1D[i] = evSys.first[i];
        weights1D[i] = 2 * std::pow( evSys.second( 0, i ), 2 );
    }

    // sort nodes increasingly and also reorder weigths for consistency
    std::vector<unsigned> sortOrder( nodes1D.size() );
    std::iota( sortOrder.begin(), sortOrder.end(), 0 );
    std::sort( sortOrder.begin(), sortOrder.end(), [&]( unsigned i, unsigned j ) { return nodes1D[i] < nodes1D[j]; } );
    Vector sorted_nodes( static_cast<unsigned>( sortOrder.size() ) ), sorted_weights( static_cast<unsigned>( sortOrder.size() ) );
    std::transform( sortOrder.begin(), sortOrder.end(), sorted_nodes.begin(), [&]( unsigned i ) { return nodes1D[i]; } );
    std::transform( sortOrder.begin(), sortOrder.end(), sorted_weights.begin(), [&]( unsigned i ) { return weights1D[i]; } );
    nodes1D   = sorted_nodes;
    weights1D = sorted_weights;

    // setup equidistant angle phi around z axis
    Vector phi( 2 * _order );
    for( unsigned i = 0; i < 2 * _order; ++i ) {
        phi[i] = ( i + 0.5 ) * M_PI / _order;
    }

    unsigned range             = std::floor( _order / 2.0 );    // comment (steffen): Only half of the points, due to projection
    double normalizationFactor = .5;

    // resize points and weights
    _pointsKarth.resize( _nq );
    _pointsSphere.resize( _nq );
    for( auto& p : _pointsKarth ) {
        p.resize( 3 );
    }
    for( auto& p : _pointsSphere ) {
        p.resize( 2 );
    }

    _weights.resize( _nq );

    // transform tensorized (x,y,z)-grid to spherical grid points
    for( unsigned j = 0; j < range; ++j ) {
        for( unsigned i = 0; i < 2 * _order; ++i ) {
            _pointsKarth[j * ( 2 * _order ) + i][0] = sqrt( 1 - nodes1D[j] * nodes1D[j] ) * std::cos( phi[i] );
            _pointsKarth[j * ( 2 * _order ) + i][1] = sqrt( 1 - nodes1D[j] * nodes1D[j] ) * std::sin( phi[i] );
            _pointsKarth[j * ( 2 * _order ) + i][2] = 0;

            _pointsSphere[j * ( 2 * _order ) + i][0] = nodes1D[j];    // my
            _pointsSphere[j * ( 2 * _order ) + i][1] = phi[i];        // phi

            _weights[j * ( 2 * _order ) + i] = normalizationFactor * M_PI / _order * weights1D[j];
        }
    }
}
