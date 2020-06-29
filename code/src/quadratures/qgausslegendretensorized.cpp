#include "../../include/quadratures/qgausslegendretensorized.h"

QGaussLegendreTensorized::QGaussLegendreTensorized( unsigned order ) : QuadratureBase( order ) {
    SetName();
    CheckOrder();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
}

void QGaussLegendreTensorized::SetPointsAndWeights() {
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
    for( unsigned i = 0; i < _order; ++i ) {
        nodes1D[i] = ( nodes1D[i] + 1.0 ) * 0.5;
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

    // only use the spheres upper half
    unsigned range = std::floor( _order / 2.0 );

    // resize points and weights
    _points.resize( _nq );
    for( auto& p : _points ) {
        p.resize( 3 );
    }
    _weights.resize( _nq );

    // transform tensorized (x,y,z)-grid to spherical grid points
    for( unsigned j = 0; j < range; ++j ) {
        for( unsigned i = 0; i < 2 * _order; ++i ) {
            _points[j * ( 2 * _order ) + i][0] = sqrt( 1 - nodes1D[j] * nodes1D[j] ) * std::cos( phi[i] );
            _points[j * ( 2 * _order ) + i][1] = sqrt( 1 - nodes1D[j] * nodes1D[j] ) * std::sin( phi[i] );
            _points[j * ( 2 * _order ) + i][2] = nodes1D[j];
            _weights[j * ( 2 * _order ) + i]   = 2.0 * M_PI / _order * weights1D[j];
        }
    }
}

void QGaussLegendreTensorized::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

std::pair<Vector, Matrix> QGaussLegendreTensorized::ComputeEigenValTriDiagMatrix( const Matrix& mat ) {
    // copied from 'Numerical Recipes' and updated + modified to work with blaze
    unsigned n = mat.rows();

    Vector d( n, 0.0 ), e( n, 0.0 );
    Matrix z( n, n, 0.0 );
    for( unsigned i = 0; i < n; ++i ) {
        d[i]          = mat( i, i );
        z( i, i )     = 1.0;
        i == 0 ? e[i] = 0.0 : e[i] = mat( i, i - 1 );
    }

    int m, l, iter, i, k;
    m = l = iter = i = k = 0;
    double s, r, p, g, f, dd, c, b;
    s = r = p = g = f = dd = c = b = 0.0;

    const double eps = std::numeric_limits<double>::epsilon();
    for( i = 1; i < static_cast<int>( n ); i++ ) e[i - 1] = e[i];
    e[n - 1] = 0.0;
    for( l = 0; l < static_cast<int>( n ); l++ ) {
        iter = 0;
        do {
            for( m = l; m < static_cast<int>( n ) - 1; m++ ) {
                dd = std::fabs( d[m] ) + std::fabs( d[m + 1] );
                if( std::fabs( e[m] ) <= eps * dd ) break;
            }
            if( m != l ) {
                if( iter++ == 30 ) ErrorMessages::Error( "Solving the tridiagonal matrix took too many iterations!", CURRENT_FUNCTION );
                g = ( d[l + 1] - d[l] ) / ( 2.0 * e[l] );
                r = Pythag( g, 1.0 );
                g = d[m] - d[l] + e[l] / ( g + std::copysign( r, g ) );
                s = c = 1.0;
                p     = 0.0;
                for( i = m - 1; i >= l; i-- ) {
                    f        = s * e[i];
                    b        = c * e[i];
                    e[i + 1] = ( r = Pythag( f, g ) );
                    if( r == 0.0 ) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s        = f / r;
                    c        = g / r;
                    g        = d[i + 1] - p;
                    r        = ( d[i] - g ) * s + 2.0 * c * b;
                    d[i + 1] = g + ( p = s * r );
                    g        = c * r - b;
                    for( k = 0; k < static_cast<int>( n ); k++ ) {
                        f = z( static_cast<unsigned>( k ), static_cast<unsigned>( i ) + 1 );
                        z( static_cast<unsigned>( k ), static_cast<unsigned>( i ) + 1 ) =
                            s * z( static_cast<unsigned>( k ), static_cast<unsigned>( i ) ) + c * f;
                        z( static_cast<unsigned>( k ), static_cast<unsigned>( i ) ) =
                            c * z( static_cast<unsigned>( k ), static_cast<unsigned>( i ) ) - s * f;
                    }
                }
                if( r == 0.0 && i >= l ) continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while( m != l );
    }
    return std::make_pair( d, z );
}

double QGaussLegendreTensorized::Pythag( const double a, const double b ) {
    // copied from 'Numerical Recipes'
    double absa = std::fabs( a ), absb = std::fabs( b );
    return ( absa > absb ? absa * std::sqrt( 1.0 + ( absb / absa ) * ( absb / absa ) )
                         : ( absb == 0.0 ? 0.0 : absb * std::sqrt( 1.0 + ( absa / absb ) * ( absa / absb ) ) ) );
}

bool QGaussLegendreTensorized::CheckOrder() {
    if( _order % 2 == 1 ) {    // order needs to be even
        ErrorMessages::Error( "ERROR! Order " + std::to_string( _order ) + " for " + GetName() + " not available. \n Order must be an even number. ",
                              CURRENT_FUNCTION );
    }
    return true;
}
