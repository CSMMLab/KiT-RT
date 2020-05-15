#include "mesh.h"

Mesh::Mesh( std::vector<Vector> nodes,
            std::vector<std::vector<unsigned>> cells,
            std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> boundaries )
    : _log( spdlog::get( "event" ) ), _dim( nodes[0].size() ), _numCells( cells.size() ), _numNodes( nodes.size() ),
      _numNodesPerCell( cells[0].size() ), _numBoundaries( boundaries.size() ), _ghostCellID( _numCells ), _nodes( nodes ), _cells( cells ),
      _boundaries( boundaries ) {
    ComputeCellAreas();
    ComputeCellMidpoints();
    ComputeConnectivity();
    ComputePartitioning();
}

Mesh::~Mesh() {}

void Mesh::ComputeConnectivity() {
    int comm_size, comm_rank;
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );
    // unsigned chunkSize    = std::ceil( static_cast<float>( _numCells ) / static_cast<float>( comm_size ) );
    // unsigned mpiCellStart = comm_rank * chunkSize;
    // unsigned mpiCellEnd   = std::min( ( comm_rank + 1 ) * chunkSize, _numCells );
    // std::vector<int> neighborsFlatPart( _numNodesPerCell * chunkSize, -1 );
    // std::vector<Vector> normalsFlatPart( _numNodesPerCell * chunkSize, Vector( _dim, -1.0 ) );

    unsigned chunkSize    = _numCells;
    unsigned mpiCellStart = 0u;
    unsigned mpiCellEnd   = _numCells;
    std::vector<int> neighborsFlatPart( _numNodesPerCell * _numCells, -1 );
    std::vector<Vector> normalsFlatPart( _numNodesPerCell * _numCells, Vector( _dim, -1.0 ) );

    auto sortedCells( _cells );
    for( unsigned i = 0; i < _numCells; ++i ) {
        std::sort( sortedCells[i].begin(), sortedCells[i].end() );
    }
    std::vector<std::vector<unsigned>> sortedBoundaries;
    for( unsigned i = 0; i < _numBoundaries; ++i ) {
        sortedBoundaries.push_back( _boundaries[i].second );
        std::sort( sortedBoundaries[i].begin(), sortedBoundaries[i].end() );
    }
    //#pragma omp parallel for
    for( unsigned i = mpiCellStart; i < mpiCellEnd; ++i ) {
        std::vector<unsigned>* cellsI = &sortedCells[i];
        for( unsigned j = 0; j < _numCells; ++j ) {
            if( i == j ) continue;
            std::vector<unsigned>* cellsJ = &sortedCells[j];
            std::vector<unsigned> commonElements;
            std::set_intersection( cellsI->begin(), cellsI->end(), cellsJ->begin(), cellsJ->end(), std::back_inserter( commonElements ) );
            if( commonElements.size() == _dim ) {
                unsigned pos0 = _numNodesPerCell * ( i - mpiCellStart );
                unsigned pos  = pos0;
                while( neighborsFlatPart[pos] != -1 && pos < pos0 + _numNodesPerCell - 1 && pos < chunkSize * _numNodesPerCell - 1 ) pos++;
                neighborsFlatPart[pos] = j;
                normalsFlatPart[pos]   = ComputeOutwardFacingNormal( _nodes[commonElements[0]], _nodes[commonElements[1]], _cellMidPoints[i] );
            }
        }
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
    std::vector<int> neighborsFlat( _numNodesPerCell * chunkSize * comm_size, -1 );
    std::vector<Vector> normalsFlat( _numNodesPerCell * chunkSize * comm_size, Vector( _dim, 0.0 ) );
    if( comm_size == 1 ) {
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

    if( std::any_of( neighborsFlat.begin(), neighborsFlat.end(), []( int i ) { return i == -1; } ) ) {
        for( unsigned idx = 0; idx < neighborsFlat.size(); ++idx ) {
            if( neighborsFlat[idx] == -1 ) _log->error( "[mesh] Detected unassigned faces at index {0}", idx );
            break;
        }
        exit( EXIT_FAILURE );
    }

    _cellNeighbors.resize( _numCells );
    _cellNormals.resize( _numCells );
    for( unsigned i = 0; i < neighborsFlat.size(); ++i ) {
        unsigned IDi = static_cast<unsigned>( i / 3.0 );
        unsigned IDj = neighborsFlat[i];
        if( IDi == IDj ) continue;
        if( IDj == _ghostCellID ) {
            if( std::find( _cellNeighbors[IDi].begin(), _cellNeighbors[IDi].end(), _ghostCellID ) == _cellNeighbors[IDi].end() ) {
                _cellNeighbors[IDi].push_back( _ghostCellID );
                _cellNormals[IDi].push_back( normalsFlat[i] );
            }
        }
        else {
            if( std::find( _cellNeighbors[IDi].begin(), _cellNeighbors[IDi].end(), IDj ) == _cellNeighbors[IDi].end() ) {
                _cellNeighbors[IDi].push_back( IDj );
                _cellNormals[IDi].push_back( normalsFlat[i] );
            }
        }
    }

    unsigned ctr = 0;
    for( unsigned i = 0; i < _numCells; ++i ) {
        if( _cellNeighbors[i].size() > _numNodesPerCell ) {
            std::cout << i << "\t\t";
            for( unsigned j = 0; j < _cellNeighbors[i].size(); ++j ) {
                std::cout << _cellNeighbors[i][j] << "\t";
            }
            std::cout << std::endl;
            ctr++;
        }
    }
    std::cout << ctr << std::endl;

    ctr = 0;
    for( unsigned i = 0; i < _numCells; ++i ) {
        if( _cellNeighbors[i].size() < _numNodesPerCell ) {
            std::cout << i << "\t\t";
            for( unsigned j = 0; j < _cellNeighbors[i].size(); ++j ) {
                std::cout << _cellNeighbors[i][j] << "\t";
            }
            std::cout << std::endl;
            ctr++;
        }
    }
    std::cout << ctr << std::endl << std::endl << std::endl;

    _isBoundaryCell.resize( _numCells, false );
    for( unsigned i = 0; i < _numCells; ++i ) {
        if( std::any_of( _cellNeighbors[i].begin(), _cellNeighbors[i].end(), [this]( unsigned i ) { return i == this->_ghostCellID; } ) )
            _isBoundaryCell[i] = true;
    }
}

void Mesh::ComputeCellAreas() {
    _cellAreas.resize( _numCells );
    for( unsigned i = 0; i < _numCells; ++i ) {
        switch( _numNodesPerCell ) {
            case 3: {
                _cellAreas[i] = std::abs( ( _nodes[_cells[i][0]][0] * ( _nodes[_cells[i][1]][1] - _nodes[_cells[i][2]][1] ) +
                                            _nodes[_cells[i][1]][0] * ( _nodes[_cells[i][2]][1] - _nodes[_cells[i][0]][1] ) +
                                            _nodes[_cells[i][2]][0] * ( _nodes[_cells[i][0]][1] - _nodes[_cells[i][1]][1] ) ) /
                                          2 );
                break;
            }
            case 4: {
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
                exit( EXIT_FAILURE );
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
        _cellMidPoints[j] = _cellMidPoints[j] / double( _cells[j].size() );
    }
}

Vector Mesh::ComputeOutwardFacingNormal( const Vector& nodeA, const Vector& nodeB, const Vector& cellCenter ) {
    double dx = nodeA[0] - nodeB[0];
    double dy = nodeA[1] - nodeB[1];
    Vector n{ -dy, dx };
    Vector p{ nodeA[0], nodeA[1] };
    if( dot( n, cellCenter - p ) > 0 ) {
        n *= -1.0;
    }
    return n;
}

void Mesh::ComputePartitioning() {
    int comm_size, comm_rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size( comm, &comm_size );
    MPI_Comm_rank( comm, &comm_rank );
    unsigned ompNThreads = omp_get_max_threads();

    if( ompNThreads > 1 ) {
        blaze::CompressedMatrix<bool> adjMatrix( _numNodes, _numNodes );
        for( unsigned i = 0; i < _numNodes; ++i ) {
            for( unsigned j = 0; j < _numCells; ++j ) {
                for( unsigned k = 0; k < _numNodesPerCell; ++k ) {
                    if( i == _cells[j][k] ) {
                        if( k == 0 ) {
                            adjMatrix.set( i, _cells[j][_numNodesPerCell - 1], true );
                            adjMatrix.set( i, _cells[j][1], true );
                        }
                        else if( k == _numNodesPerCell - 1 ) {
                            adjMatrix.set( i, _cells[j][_numNodesPerCell - 2], true );
                            adjMatrix.set( i, _cells[j][0], true );
                        }
                        else {
                            adjMatrix.set( i, _cells[j][k - 1], true );
                            adjMatrix.set( i, _cells[j][k + 1], true );
                        }
                    }
                }
            }
        }

        std::vector<int> xadj( _numNodes + 1 );
        int ctr = 0;
        std::vector<int> adjncy;
        for( unsigned i = 0; i < _numNodes; ++i ) {
            xadj[i] = ctr;
            for( unsigned j = 0; j < _numNodes; ++j ) {
                if( adjMatrix( i, j ) ) {
                    adjncy.push_back( static_cast<int>( j ) );
                    ctr++;
                }
            }
        }
        xadj[_numNodes] = ctr;

        std::vector<int> partitions;
        int edgecut  = 0;
        int ncon     = 1;
        int nparts   = ompNThreads;
        real_t ubvec = static_cast<real_t>( 1.05 );

        if( comm_size > 1 ) {    // use parmetis
            std::vector<int> local_xadj{ 0 };
            unsigned xadjChunk = std::ceil( static_cast<float>( _numNodes ) / static_cast<float>( comm_size ) );
            unsigned xadjStart = comm_rank * xadjChunk + 1;
            unsigned adjncyStart;
            for( unsigned i = 0; i < xadjChunk; ++i ) {
                if( i == 0 ) adjncyStart = xadj[i + xadjStart - 1];
                local_xadj.push_back( xadj[i + xadjStart] );
            }
            unsigned adjncyEnd = local_xadj.back();
            for( unsigned i = 1; i < local_xadj.size(); ++i ) {
                local_xadj[i] -= xadj[xadjStart - 1];
            }
            std::vector<int> local_adjncy( adjncy.begin() + adjncyStart, adjncy.begin() + adjncyEnd );

            std::vector<int> vtxdist{ 0 };
            for( unsigned i = 1; i <= static_cast<unsigned>( comm_size ); i++ ) {
                vtxdist.push_back( std::min( i * xadjChunk, _numNodes ) );
            }

            real_t* tpwgts = new real_t[nparts];
            for( unsigned i = 0; i < static_cast<unsigned>( nparts ); ++i ) {
                tpwgts[i] = 1.0 / static_cast<real_t>( nparts );
            }

            int wgtflag = 0;
            int numflag = 0;

            int options[1];
            options[0] = 0;

            unsigned chunkSize = local_xadj.size() - 1;
            int* part          = new int[chunkSize];
            ParMETIS_V3_PartKway( vtxdist.data(),
                                  local_xadj.data(),
                                  local_adjncy.data(),
                                  nullptr,
                                  nullptr,
                                  &wgtflag,
                                  &numflag,
                                  &ncon,
                                  &nparts,
                                  tpwgts,
                                  &ubvec,
                                  options,
                                  &edgecut,
                                  part,
                                  &comm );
            partitions.resize( chunkSize * comm_size );
            MPI_Allgather( part, chunkSize, MPI_INT, partitions.data(), chunkSize, MPI_INT, comm );
            partitions.resize( _numNodes );

            delete[] tpwgts;
        }
        else {    // use metis
            int options[METIS_NOPTIONS];
            METIS_SetDefaultOptions( options );
            int nvtxs = _numNodes;
            int vsize;
            partitions.resize( _numNodes );
            METIS_PartGraphKway(
                &nvtxs, &ncon, xadj.data(), adjncy.data(), nullptr, &vsize, nullptr, &nparts, nullptr, &ubvec, options, &edgecut, partitions.data() );
        }
        _colors.resize( _numCells );
        for( unsigned i = 0; i < _numCells; ++i ) {
            std::map<unsigned, int> occurances;
            for( unsigned j = 0; j < _numNodesPerCell; ++j ) {
                ++occurances[partitions[_cells[i][j]]];
            }
            _colors[i] = occurances.rbegin()->first;
        }
    }
    else {
        _colors.resize( _numCells, 0u );
    }
}

void Mesh::ComputeSlopes( unsigned nq, VectorVector& psiDerX, VectorVector& psiDerY, const VectorVector& psi ) const {
    for( unsigned k = 0; k < nq; ++k ) {
        for( unsigned j = 0; j < _numCells; ++j ) {
            // if( cell->IsBoundaryCell() ) continue; // skip ghost cells
            // compute derivative by summing over cell boundary
            for( unsigned l = 0; l < _cellNeighbors[j].size(); ++l ) {
                psiDerX[j][k] = psiDerX[j][k] + 0.5 * ( psi[j][k] + psi[_cellNeighbors[j][l]][k] ) * _cellNormals[j][l][0] / _cellAreas[j];
                psiDerY[j][k] = psiDerY[j][k] + 0.5 * ( psi[j][k] + psi[_cellNeighbors[j][l]][k] ) * _cellNormals[j][l][1] / _cellAreas[j];
            }
        }
    }
}

unsigned Mesh::GetDim() const { return _dim; }
unsigned Mesh::GetNumCells() const { return _numCells; }
unsigned Mesh::GetNumNodes() const { return _numNodes; }
unsigned Mesh::GetNumNodesPerCell() const { return _numNodesPerCell; }

const std::vector<Vector>& Mesh::GetNodes() const { return _nodes; }
const std::vector<Vector>& Mesh::GetCellMidPoints() const { return _cellMidPoints; }
const std::vector<std::vector<unsigned>>& Mesh::GetCells() const { return _cells; }
const std::vector<double>& Mesh::GetCellAreas() const { return _cellAreas; }
const std::vector<unsigned>& Mesh::GetPartitionIDs() const { return _colors; }
const std::vector<std::vector<unsigned>>& Mesh::GetNeighbours() const { return _cellNeighbors; }
const std::vector<std::vector<Vector>>& Mesh::GetNormals() const { return _cellNormals; }
const std::vector<bool>& Mesh::GetBoundaryCellArray() const { return _isBoundaryCell; }
