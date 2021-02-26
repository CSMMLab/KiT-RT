#include "common/mesh.h"

#include <chrono>

Mesh::Mesh( std::vector<Vector> nodes,
            std::vector<std::vector<unsigned>> cells,
            std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> boundaries )
    : _dim( nodes[0].size() ), _numCells( cells.size() ), _numNodes( nodes.size() ), _numNodesPerCell( cells[0].size() ),
      _numBoundaries( boundaries.size() ), _ghostCellID( _numCells ), _nodes( nodes ), _cells( cells ), _boundaries( boundaries ) {
    if( _dim == 2 ) {
        _numNodesPerBoundary = 2u;
    }
    else {
        ErrorMessages::Error( "Unsupported mesh dimension!", CURRENT_FUNCTION );
    }

    ComputeCellAreas();
    ComputeCellMidpoints();
    ComputeConnectivity();
    ComputeBounds();
}

Mesh::~Mesh() {}

void Mesh::ComputeConnectivity() {

    // get MPI info
    int comm_size, comm_rank;
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );

    // determine number/chunk size and indices of cells treated by each mpi thread
    unsigned chunkSize    = std::ceil( static_cast<float>( _numCells ) / static_cast<float>( comm_size ) );
    unsigned mpiCellStart = comm_rank * chunkSize;
    unsigned mpiCellEnd   = std::min( ( comm_rank + 1 ) * chunkSize, _numCells );

    // 'flat' vectors are a flattened representation of the neighbors numCells<numNodesPerCell> nested vectors; easier for MPI
    // 'part' vectors store information for each single MPI thread
    std::vector<int> neighborsFlatPart( _numNodesPerCell * chunkSize, -1 );
    std::vector<Vector> normalsFlatPart( _numNodesPerCell * chunkSize, Vector( _dim, -1.0 ) );

    // pre sort cells and boundaries; sorting is needed for std::set_intersection
    auto sortedCells( _cells );
    for( unsigned i = 0; i < _numCells; ++i ) {
        std::sort( sortedCells[i].begin(), sortedCells[i].end() );
    }
    std::vector<std::vector<unsigned>> sortedBoundaries;
    for( unsigned i = 0; i < _numBoundaries; ++i ) {
        sortedBoundaries.push_back( _boundaries[i].second );
        std::sort( sortedBoundaries[i].begin(), sortedBoundaries[i].end() );
    }

    blaze::CompressedMatrix<bool> connMat( _numCells, _numNodes );
    for( unsigned i = mpiCellStart; i < mpiCellEnd; ++i ) {
        for( auto j : _cells[i] ) connMat.set( i, j, true );
    }

// determine neighbor cells and normals with MPI and OpenMP
#pragma omp parallel for
    for( unsigned i = mpiCellStart; i < mpiCellEnd; ++i ) {
        std::vector<unsigned>* cellsI = &sortedCells[i];
        for( unsigned j = 0; j < _numCells; ++j ) {
            if( i == j ) continue;
            if( static_cast<unsigned>( blaze::dot( blaze::row( connMat, i ), blaze::row( connMat, j ) ) ) ==
                _numNodesPerBoundary ) {    // in 2D cells are neighbors if they share two nodes std::vector<unsigned>* cellsJ = &sortedCells[j];
                std::vector<unsigned>* cellsJ = &sortedCells[j];
                std::vector<unsigned> commonElements;
                std::set_intersection( cellsI->begin(),
                                       cellsI->end(),
                                       cellsJ->begin(),
                                       cellsJ->end(),
                                       std::back_inserter( commonElements ) );    // find common nodes of two cells
                // determine unused index
                unsigned pos0 = _numNodesPerCell * ( i - mpiCellStart );
                unsigned pos  = pos0;
                while( neighborsFlatPart[pos] != -1 && pos < pos0 + _numNodesPerCell - 1 && pos < chunkSize * _numNodesPerCell - 1 ) pos++;
                neighborsFlatPart[pos] = j;
                // compute normal vector
                normalsFlatPart[pos] = ComputeOutwardFacingNormal( _nodes[commonElements[0]], _nodes[commonElements[1]], _cellMidPoints[i] );
            }
        }
        // boundaries are treated similarly to normal cells, but need a special treatment due to the absence of a neighboring cell
        for( unsigned k = 0; k < _boundaries.size(); ++k ) {
            std::vector<unsigned>* bNodes = &sortedBoundaries[k];
            for( unsigned j = 0; j < _boundaries[k].second.size(); ++j ) {
                std::vector<unsigned> commonElements;
                std::set_intersection( cellsI->begin(), cellsI->end(), bNodes->begin(), bNodes->end(), std::back_inserter( commonElements ) );
                if( commonElements.size() == _dim ) {
                    unsigned pos0 = _numNodesPerCell * ( i - mpiCellStart );
                    unsigned pos  = pos0;
                    while( neighborsFlatPart[pos] != -1 && pos < pos0 + _numNodesPerCell - 1 && pos < chunkSize * _numNodesPerCell - 1 ) pos++;
                    neighborsFlatPart[pos] = _ghostCellID;
                    normalsFlatPart[pos]   = ComputeOutwardFacingNormal( _nodes[commonElements[0]], _nodes[commonElements[1]], _cellMidPoints[i] );
                }
            }
        }
    }

    // gather distributed data on all MPI threads
    std::vector<int> neighborsFlat( _numNodesPerCell * chunkSize * comm_size, -1 );
    std::vector<Vector> normalsFlat( _numNodesPerCell * chunkSize * comm_size, Vector( _dim, 0.0 ) );
    if( comm_size == 1 ) {    // can be done directly if there is only one MPI thread
        neighborsFlat.assign( neighborsFlatPart.begin(), neighborsFlatPart.end() );
        normalsFlat.assign( normalsFlatPart.begin(), normalsFlatPart.end() );
    }
    else {
        MPI_Allgather( neighborsFlatPart.data(),
                       _numNodesPerCell * chunkSize,
                       MPI_INT,
                       neighborsFlat.data(),
                       _numNodesPerCell * chunkSize,
                       MPI_INT,
                       MPI_COMM_WORLD );
    }

    // check for any unassigned faces
    if( std::any_of( neighborsFlat.begin(), neighborsFlat.end(), []( int i ) { return i == -1; } ) ) {
        for( unsigned idx = 0; idx < neighborsFlat.size(); ++idx ) {
            if( neighborsFlat[idx] == -1 )
                ErrorMessages::Error( "Detected unassigned faces at index " + std::to_string( idx ) + " !", CURRENT_FUNCTION );
        }
    }

    // reorder neighbors and normals into nested structure
    _cellNeighbors.resize( _numCells );
    _cellNormals.resize( _numCells );
    for( unsigned i = 0; i < neighborsFlat.size(); ++i ) {
        unsigned IDi = static_cast<unsigned>( i / static_cast<double>( _numNodesPerCell ) );
        unsigned IDj = neighborsFlat[i];
        if( IDi == IDj ) continue;     // avoid self assignment
        if( IDj == _ghostCellID ) {    // cell is boundary cell
            if( std::find( _cellNeighbors[IDi].begin(), _cellNeighbors[IDi].end(), _ghostCellID ) == _cellNeighbors[IDi].end() ) {
                _cellNeighbors[IDi].push_back( _ghostCellID );
                _cellNormals[IDi].push_back( normalsFlat[i] );
            }
        }
        else {    // normal cell neighbor
            if( std::find( _cellNeighbors[IDi].begin(), _cellNeighbors[IDi].end(), IDj ) == _cellNeighbors[IDi].end() ) {
                _cellNeighbors[IDi].push_back( IDj );
                _cellNormals[IDi].push_back( normalsFlat[i] );
            }
        }
    }

    // assign boundary types to all cells
    _cellBoundaryTypes.resize( _numCells, BOUNDARY_TYPE::NONE );
    for( unsigned i = 0; i < _numCells; ++i ) {
        if( std::any_of( _cellNeighbors[i].begin(), _cellNeighbors[i].end(), [this]( unsigned i ) {
                return i == _ghostCellID;
            } ) ) {                           // cell is boundary cell
            for( auto bc : _boundaries ) {    // loop over all boundaries and boundary nodes and search for the respective nodeID
                for( auto cNodes : _cells[i] ) {
                    if( std::find( bc.second.begin(), bc.second.end(), cNodes ) != bc.second.end() ) {
                        _cellBoundaryTypes[i] = bc.first;
                    }
                }
            }
        }
    }
}

void Mesh::ComputeCellAreas() {
    _cellAreas.resize( _numCells );
    for( unsigned i = 0; i < _numCells; ++i ) {
        switch( _numNodesPerCell ) {
            case 3: {    // triangular cells
                _cellAreas[i] = std::abs( ( _nodes[_cells[i][0]][0] * ( _nodes[_cells[i][1]][1] - _nodes[_cells[i][2]][1] ) +
                                            _nodes[_cells[i][1]][0] * ( _nodes[_cells[i][2]][1] - _nodes[_cells[i][0]][1] ) +
                                            _nodes[_cells[i][2]][0] * ( _nodes[_cells[i][0]][1] - _nodes[_cells[i][1]][1] ) ) /
                                          2 );
                break;
            }
            case 4: {    // quadrilateral cells
                std::vector d1{ _nodes[_cells[i][0]][0] - _nodes[_cells[i][1]][0], _nodes[_cells[i][0]][1] - _nodes[_cells[i][1]][1] };
                std::vector d2{ _nodes[_cells[i][1]][0] - _nodes[_cells[i][2]][0], _nodes[_cells[i][1]][1] - _nodes[_cells[i][2]][1] };
                std::vector d3{ _nodes[_cells[i][2]][0] - _nodes[_cells[i][3]][0], _nodes[_cells[i][2]][1] - _nodes[_cells[i][3]][1] };
                std::vector d4{ _nodes[_cells[i][3]][0] - _nodes[_cells[i][0]][0], _nodes[_cells[i][3]][1] - _nodes[_cells[i][0]][1] };

                double a = std::sqrt( d1[0] * d1[0] + d1[1] * d1[1] );
                double b = std::sqrt( d2[0] * d2[0] + d2[1] * d2[1] );
                double c = std::sqrt( d3[0] * d3[0] + d3[1] * d3[1] );
                double d = std::sqrt( d4[0] * d4[0] + d4[1] * d4[1] );
                double T = 0.5 * ( a + b + c + d );

                double alpha = std::acos( ( d4[0] * d1[0] + d4[1] * d1[1] ) / ( a * d ) );
                double beta  = std::acos( ( d2[0] * d3[0] + d2[1] * d3[1] ) / ( b * c ) );

                _cellAreas[i] = std::sqrt( ( T - a ) * ( T - b ) * ( T - c ) * ( T - d ) -
                                           a * b * c * d * std::cos( 0.5 * ( alpha + beta ) ) * std::cos( 0.5 * ( alpha + beta ) ) );
                break;
            }
            default: {
                ErrorMessages::Error( "Area computation for cells with " + std::to_string( _numNodesPerCell ) + " nodes is not implemented yet!",
                                      CURRENT_FUNCTION );
            }
        }
    }
}

void Mesh::ComputeCellMidpoints() {
    _cellMidPoints = std::vector( _numCells, Vector( _dim, 0.0 ) );
    for( unsigned j = 0; j < _numCells; ++j ) {
        for( unsigned l = 0; l < _cells[j].size(); ++l ) {
            _cellMidPoints[j] = _cellMidPoints[j] + _nodes[_cells[j][l]];
        }
        _cellMidPoints[j] = _cellMidPoints[j] / static_cast<double>( _cells[j].size() );
    }
}

Vector Mesh::ComputeOutwardFacingNormal( const Vector& nodeA, const Vector& nodeB, const Vector& cellCenter ) {
    double dx = nodeA[0] - nodeB[0];
    double dy = nodeA[1] - nodeB[1];
    Vector n{ -dy, dx };                    // normal vector
    Vector p{ nodeA[0], nodeA[1] };         // take arbitrary node
    if( dot( n, cellCenter - p ) > 0 ) {    // dot product is negative for outward facing normals
        n *= -1.0;                          // if dot product is positive -> flip normal
    }
    return n;
}

void Mesh::ComputeSlopes( unsigned nq, VectorVector& psiDerX, VectorVector& psiDerY, const VectorVector& psi ) const {
    for( unsigned k = 0; k < nq; ++k ) {
        for( unsigned j = 0; j < _numCells; ++j ) {
            psiDerX[j][k] = 0.0;
            psiDerY[j][k] = 0.0;

            // if( cell->IsBoundaryCell() ) continue; // skip ghost cells
            if( _cellBoundaryTypes[j] != 2 ) continue;    // skip ghost cells
            // compute derivative by summing over cell boundary
            for( unsigned l = 0; l < _cellNeighbors[j].size(); ++l ) {
                psiDerX[j][k] = psiDerX[j][k] + 0.5 * ( psi[j][k] + psi[_cellNeighbors[j][l]][k] ) * _cellNormals[j][l][0] / _cellAreas[j];
                psiDerY[j][k] = psiDerY[j][k] + 0.5 * ( psi[j][k] + psi[_cellNeighbors[j][l]][k] ) * _cellNormals[j][l][1] / _cellAreas[j];
            }
        }
    }
}

void Mesh::ReconstructSlopesS( unsigned nq, VectorVector& psiDerX, VectorVector& psiDerY, const VectorVector& psi ) const {

    double duxL;
    double duyL;
    double duxR;
    double duyR;

    for( unsigned k = 0; k < nq; ++k ) {
        for( unsigned j = 0; j < _numCells; ++j ) {
            // reset derivatives
            psiDerX[j][k] = 0.0;
            psiDerY[j][k] = 0.0;

            // skip boundary cells
            if( _cellBoundaryTypes[j] != 2 ) continue;

            duxL = ( psi[j][k] - psi[_cellNeighbors[j][0]][k] ) / ( _cellMidPoints[j][0] - _cellMidPoints[_cellNeighbors[j][0]][0] + 1e-8 );
            duyL = ( psi[j][k] - psi[_cellNeighbors[j][0]][k] ) / ( _cellMidPoints[j][1] - _cellMidPoints[_cellNeighbors[j][0]][1] + 1e-8 );

            duxR = ( psi[j][k] - psi[_cellNeighbors[j][2]][k] ) / ( _cellMidPoints[j][0] - _cellMidPoints[_cellNeighbors[j][2]][0] + 1e-8 );
            duyR = ( psi[j][k] - psi[_cellNeighbors[j][2]][k] ) / ( _cellMidPoints[j][1] - _cellMidPoints[_cellNeighbors[j][2]][1] + 1e-8 );

            psiDerX[j][k] = LMinMod( duxL, duxR ) * 0.5;
            psiDerY[j][k] = LMinMod( duyL, duyR ) * 0.5;
        }
    }
}

void Mesh::ReconstructSlopesU( unsigned nq, VectorVector& psiDerX, VectorVector& psiDerY, const VectorVector& psi ) const {

    double phi;
    // VectorVector dPsiMax = std::vector( _numCells, Vector( nq, 0.0 ) );
    // VectorVector dPsiMin = std::vector( _numCells, Vector( nq, 0.0 ) );
    VectorVector dPsiMax( _numCells, Vector( nq, 0.0 ) );
    VectorVector dPsiMin( _numCells, Vector( nq, 0.0 ) );

    std::vector<std::vector<Vector>> psiSample( _numCells, std::vector<Vector>( nq, Vector( _numNodesPerCell, 0.0 ) ) );
    std::vector<std::vector<Vector>> phiSample( _numCells, std::vector<Vector>( nq, Vector( _numNodesPerCell, 0.0 ) ) );

    for( unsigned k = 0; k < nq; ++k ) {
        for( unsigned j = 0; j < _numCells; ++j ) {
            // reset derivatives
            psiDerX[j][k] = 0.0;
            psiDerY[j][k] = 0.0;

            // skip boundary cells
            if( _cellBoundaryTypes[j] != 2 ) continue;

            /*
            // version 1: original taste
            // step 1: calculate psi difference around neighbors and theoretical derivatives by Gauss theorem
            for( unsigned l = 0; l < _cellNeighbors[j].size(); ++l ) {
                if( psi[_cellNeighbors[j][l]][k] - psi[j][k] > dPsiMax[j][k] ) dPsiMax[j][k] = psi[_cellNeighbors[j][l]][k] - psi[j][k];
                if( psi[_cellNeighbors[j][l]][k] - psi[j][k] < dPsiMin[j][k] ) dPsiMin[j][k] = psi[_cellNeighbors[j][l]][k] - psi[j][k];

                psiDerX[j][k] += 0.5 * ( psi[j][k] + psi[_cellNeighbors[j][l]][k] ) * _cellNormals[j][l][0] / _cellAreas[j];
                psiDerY[j][k] += 0.5 * ( psi[j][k] + psi[_cellNeighbors[j][l]][k] ) * _cellNormals[j][l][1] / _cellAreas[j];
            }

            for( unsigned l = 0; l < _cellNeighbors[j].size(); ++l ) {
                // step 2: choose sample points
                // psiSample[j][k][l] = 0.5 * ( psi[j][k] + psi[_cellNeighbors[j][l]][k] ); // interface central points
                psiSample[j][k][l] = psi[j][k] +
                                     psiDerX[j][k] * ( _nodes[_cells[j][l]][0] - _cellMidPoints[j][0] ) +
                                     psiDerY[j][k] * ( _nodes[_cells[j][l]][1] - _cellMidPoints[j][1] );    // vertex points

                // step 3: calculate Phi_ij at sample points
                if( psiSample[j][k][l] > psi[j][k] ) {
                    phiSample[j][k][l] = fmin( 1.0, dPsiMax[j][k] / ( psiSample[j][k][l] - psi[j][k] ) );
                }
                else if( psiSample[j][l][k] < psi[j][k] ) {
                    phiSample[j][k][l] = fmin( 1.0, dPsiMin[j][k] / ( psiSample[j][k][l] - psi[j][k] ) );
                }
                else {
                    phiSample[j][k][l] = 1.0;
                }
            }

            // step 4: find minimum limiter function phi
            phi = min( phiSample[j][k] );

            // step 5: limit the slope reconstructed from Gauss theorem
            psiDerX[j][k] *= phi;
            psiDerY[j][k] *= phi;
            */

            double eps = 1e-6;
            // version 2: Venkatakrishnan limiter
            // step 1: calculate psi difference around neighbors and theoretical derivatives by Gauss theorem
            for( unsigned l = 0; l < _cellNeighbors[j].size(); ++l ) {
                if( psi[_cellNeighbors[j][l]][k] - psi[j][k] > dPsiMax[j][k] ) dPsiMax[j][k] = psi[_cellNeighbors[j][l]][k] - psi[j][k];
                if( psi[_cellNeighbors[j][l]][k] - psi[j][k] < dPsiMin[j][k] ) dPsiMin[j][k] = psi[_cellNeighbors[j][l]][k] - psi[j][k];

                psiDerX[j][k] += 0.5 * ( psi[j][k] + psi[_cellNeighbors[j][l]][k] ) * _cellNormals[j][l][0] / _cellAreas[j];
                psiDerY[j][k] += 0.5 * ( psi[j][k] + psi[_cellNeighbors[j][l]][k] ) * _cellNormals[j][l][1] / _cellAreas[j];
            }

            for( unsigned l = 0; l < _cellNeighbors[j].size(); ++l ) {
                // step 2: choose sample points
                psiSample[j][k][l] = 10.0 * psiDerX[j][k] * ( _cellMidPoints[_cellNeighbors[j][l]][0] - _cellMidPoints[j][0] ) +
                                     10.0 * psiDerY[j][k] * ( _cellMidPoints[_cellNeighbors[j][l]][1] - _cellMidPoints[j][1] );

                // step 3: calculate Phi_ij at sample points
                if( psiSample[j][k][l] > 0.0 ) {
                    phiSample[j][k][l] =
                        ( dPsiMax[j][k] * dPsiMax[j][k] + 2.0 * dPsiMax[j][k] * psiSample[j][k][l] + eps ) /
                        ( dPsiMax[j][k] * dPsiMax[j][k] + dPsiMax[j][k] * psiSample[j][k][l] + 2.0 * psiSample[j][k][l] * psiSample[j][k][l] + eps );
                }
                else if( psiSample[j][k][l] < 0.0 ) {
                    phiSample[j][k][l] =
                        ( dPsiMin[j][k] * dPsiMin[j][k] + 2.0 * dPsiMin[j][k] * psiSample[j][k][l] + eps ) /
                        ( dPsiMin[j][k] * dPsiMin[j][k] + dPsiMin[j][k] * psiSample[j][k][l] + 2.0 * psiSample[j][k][l] * psiSample[j][k][l] + eps );
                    ;
                }
                else {
                    phiSample[j][k][l] = 1.0;
                }
            }

            // step 4: find minimum limiter function phi
            phi = min( phiSample[j][k] );
            // phi = fmin( 0.5, abs(phi) );

            // step 5: limit the slope reconstructed from Gauss theorem
            psiDerX[j][k] *= phi;
            psiDerY[j][k] *= phi;
        }
    }
}

void Mesh::ComputeBounds() {
    _bounds = std::vector( _dim, std::make_pair( std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() ) );
    for( unsigned i = 0; i < _numNodes; ++i ) {
        for( unsigned j = 0; j < _dim; ++j ) {
            if( _nodes[i][j] < _bounds[j].first ) _bounds[j].first = _nodes[i][j];
            if( _nodes[i][j] > _bounds[j].second ) _bounds[j].second = _nodes[i][j];
        }
    }
}

const std::vector<Vector>& Mesh::GetNodes() const { return _nodes; }
const std::vector<Vector>& Mesh::GetCellMidPoints() const { return _cellMidPoints; }
const std::vector<std::vector<unsigned>>& Mesh::GetCells() const { return _cells; }
const std::vector<double>& Mesh::GetCellAreas() const { return _cellAreas; }
const std::vector<std::vector<unsigned>>& Mesh::GetNeighbours() const { return _cellNeighbors; }
const std::vector<std::vector<Vector>>& Mesh::GetNormals() const { return _cellNormals; }
const std::vector<BOUNDARY_TYPE>& Mesh::GetBoundaryTypes() const { return _cellBoundaryTypes; }
const std::vector<std::pair<double, double>> Mesh::GetBounds() const { return _bounds; }

double Mesh::GetDistanceToOrigin( unsigned idx_cell ) const {
    double distance = 0.0;
    for( unsigned idx_dim = 0; idx_dim < _dim; idx_dim++ ) {
        distance += _cellMidPoints[idx_cell][idx_dim] * _cellMidPoints[idx_cell][idx_dim];
    }
    return sqrt( distance );
}
