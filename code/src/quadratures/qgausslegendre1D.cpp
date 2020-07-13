#include "quadratures/qgausslegendre1D.h"
#include "toolboxes/errormessages.h"

QGaussLegendre1D::QGaussLegendre1D( unsigned order ) : QuadratureBase( order ) {
    SetName();
    CheckOrder();
    SetNq();
    SetPointsAndWeights();
    SetConnectivity();
}

void QGaussLegendre1D::SetPointsAndWeights() {
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

    // resize points and weights
    _points.resize( _nq );
    _weights.resize( _nq );
    unsigned dim = 3;
    for( unsigned k = 0; k < _nq; ++k ) {
        _points[k].resize( dim );
        _points[k][0] = nodes1D[k];
        _weights[k]   = weights1D[k];
    }
}

void QGaussLegendre1D::SetConnectivity() {    // TODO
    // Not initialized for this quadrature.
    VectorVectorU connectivity;
    _connectivity = connectivity;
}

std::pair<Vector, Matrix> QGaussLegendre1D::ComputeEigenValTriDiagMatrix( const Matrix& mat ) {
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

double QGaussLegendre1D::Pythag( const double a, const double b ) {
    // copied from 'Numerical Recipes'
    double absa = std::fabs( a ), absb = std::fabs( b );
    return ( absa > absb ? absa * std::sqrt( 1.0 + ( absb / absa ) * ( absb / absa ) )
                         : ( absb == 0.0 ? 0.0 : absb * std::sqrt( 1.0 + ( absa / absb ) * ( absa / absb ) ) ) );
}

bool QGaussLegendre1D::CheckOrder() { return true; }
