#include "mesh.h"

Mesh::Mesh( std::vector<Vector> nodes,
            std::vector<std::vector<unsigned>> cells,
            std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> boundaries )
    : _log( spdlog::get( "event" ) ), _dim( nodes[0].size() ), _numCells( cells.size() ), _numNodes( nodes.size() ),
      _numNodesPerCell( cells[0].size() ), _numBoundaries( boundaries.size() ), _ghostCellID( _numCells ), _nodes( nodes ), _cells( cells ),
      _boundaries( boundaries ) {
    ComputeCellAreas();
    ComputeNormals();
    ComputeConnectivity();
    ComputePartitioning();
}

Mesh::~Mesh() {}

void Mesh::ComputeConnectivity() {
    int comm_size, comm_rank;
    MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );
    unsigned chunkSize    = std::ceil( static_cast<float>( _numCells ) / static_cast<float>( comm_size ) );
    unsigned mpiCellStart = comm_rank * chunkSize;
    unsigned mpiCellEnd   = std::max( comm_rank + 1 * chunkSize, _numCells - 1 );
    std::vector<int> neighborsFlatPart( _numNodesPerCell * chunkSize, -1 );

    auto sortedCells( _cells );
    for( unsigned i = 0; i < _numCells; ++i ) {
        std::sort( sortedCells[i].begin(), sortedCells[i].end() );
    }
    std::vector<std::vector<unsigned>> sortedBoundaries;
    for( unsigned i = 0; i < _numBoundaries; ++i ) {
        sortedBoundaries.push_back( _boundaries[i].second );
        std::sort( sortedBoundaries[i].begin(), sortedBoundaries[i].end() );
    }
#pragma omp parallel for
    for( unsigned i = mpiCellStart; i < mpiCellEnd; ++i ) {
        std::vector<unsigned>* cellsI = &sortedCells[i];
        for( unsigned j = 0; j < _numCells; ++j ) {
            std::vector<unsigned>* cellsJ = &sortedCells[j];
            std::vector<unsigned> commonElements;
            std::set_intersection( cellsI->begin(), cellsI->end(), cellsJ->begin(), cellsJ->end(), std::back_inserter( commonElements ) );
            if( commonElements.size() == _dim ) {
                unsigned pos0 = _numNodesPerCell * ( i - mpiCellStart );
                unsigned pos  = pos0;
                while( neighborsFlatPart[pos] != -1 && pos < pos0 + _numNodesPerCell && pos < _numCells * _numNodesPerCell - 1 ) pos++;
                neighborsFlatPart[pos] = j;
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
                    while( neighborsFlatPart[pos] != -1 && pos < pos0 + _numNodesPerCell && pos < _numCells * _numNodesPerCell - 1 ) pos++;
                    neighborsFlatPart[pos] = _ghostCellID;
                }
            }
        }
    }
    std::vector<int> neighborsFlat( _numNodesPerCell * _numCells, -1 );
    if( comm_size == 1 )
        neighborsFlat.assign( neighborsFlatPart.begin(), neighborsFlatPart.end() );
    else
        MPI_Allgather(
            neighborsFlatPart.data(), chunkSize, MPI::UNSIGNED, neighborsFlat.data(), _numNodesPerCell * _numCells, MPI::UNSIGNED, MPI_COMM_WORLD );
    if( std::any_of( neighborsFlat.begin(), neighborsFlat.end(), []( int i ) { return i == -1; } ) ) {
        unsigned idx = 0;
        for( ; idx < neighborsFlat.size(); ++idx ) {
            if( neighborsFlat[idx] == -1 ) _log->error( "[mesh] Detected unassigned faces at index {0}", idx );
        }
        exit( EXIT_FAILURE );
    }
    for( auto it = neighborsFlat.begin(); it != neighborsFlat.end(); it += _numNodesPerCell )
        _cellNeighbors.push_back( std::vector<unsigned>( it, it + _numNodesPerCell ) );
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
                std::vector d1{_nodes[_cells[i][0]][0] - _nodes[_cells[i][1]][0], _nodes[_cells[i][0]][1] - _nodes[_cells[i][1]][1]};
                std::vector d2{_nodes[_cells[i][1]][0] - _nodes[_cells[i][2]][0], _nodes[_cells[i][1]][1] - _nodes[_cells[i][2]][1]};
                std::vector d3{_nodes[_cells[i][2]][0] - _nodes[_cells[i][3]][0], _nodes[_cells[i][2]][1] - _nodes[_cells[i][3]][1]};
                std::vector d4{_nodes[_cells[i][3]][0] - _nodes[_cells[i][0]][0], _nodes[_cells[i][3]][1] - _nodes[_cells[i][0]][1]};

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

Vector Mesh::ComputeOutwardFacingNormal( const Vector& nodeA, const Vector& nodeB, const Vector& cellCenter ) {
    double dx = nodeA[0] - nodeB[0];
    double dy = nodeA[1] - nodeB[1];
    Vector n{-dy, dx};
    Vector p{nodeA[0], nodeA[1]};
    if( dot( n, cellCenter - p ) > 0 ) {
        n *= -1.0;
    }
    return n;
}

void Mesh::ComputeNormals() {
    _cellNormals.resize( _numCells, std::vector<Vector>( _numNodesPerCell, Vector( _dim, 0.0 ) ) );
    for( unsigned i = 0; i < _numCells; ++i ) {
        auto cellNodes = _cells[i];
        Vector cellCenter( _dim, 0.0 );
        for( unsigned j = 0; j < _numNodesPerCell; ++j ) {
            for( unsigned k = 0; k < _dim; ++k ) {
                cellCenter[k] += _nodes[cellNodes[j]][k];
            }
        }
        for( unsigned j = 0; j < _numNodesPerCell - 1; ++j ) {
            _cellNormals[i][j] = ComputeOutwardFacingNormal( _nodes[cellNodes[j]], _nodes[cellNodes[j + 1]], cellCenter );
        }
        _cellNormals[i][_numNodesPerCell - 1] =
            ComputeOutwardFacingNormal( _nodes[cellNodes[_numNodesPerCell - 1]], _nodes[cellNodes[0]], cellCenter );
    }
}

void Mesh::ComputePartitioning() {
    int comm_size, comm_rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size( comm, &comm_size );
    MPI_Comm_rank( comm, &comm_rank );
    unsigned ompNThreads = omp_get_max_threads();

    unsigned long nPoint = _numNodes;

    if( ompNThreads > 1 ) {
        blaze::CompressedMatrix<unsigned> adjMatrix( _numNodes, _numNodes );

        int* xadj = new int[_nodes.size() + 1];
        std::vector<int> tmp_adjncy;
        int ctr = 0;
        for( unsigned i = 0; i < _numNodes; ++i ) {
            xadj[i] = ctr;
            for( unsigned j = 0; j < _numNodes; ++j ) {
                if( adjMatrix( i, j ) == 1u ) {
                    tmp_adjncy.push_back( static_cast<int>( j ) );
                    ctr++;
                }
            }
        }
        int* adjncy = new int[static_cast<unsigned long>( xadj[_numNodes + 1] )];
        std::copy( tmp_adjncy.begin(), tmp_adjncy.end(), adjncy );

        if( comm_size > 1 ) {    // use parmetis
            std::vector<unsigned> npoint_procs( ompNThreads );
            for( unsigned i = 0; i < npoint_procs.size() - 1u; ++i )
                npoint_procs[i] = static_cast<unsigned>( std::round( static_cast<double>( _nodes.size() ) / static_cast<double>( ompNThreads ) ) );
            npoint_procs.back() = static_cast<unsigned>( _numNodes ) - ( ompNThreads - 1u ) * npoint_procs[0];

            std::vector<unsigned> starting_node( ompNThreads, 0u );
            std::vector<unsigned> ending_node( ompNThreads, 0u );
            nPoint           = npoint_procs[comm_rank];
            starting_node[0] = 0;
            ending_node[0]   = starting_node[0] + npoint_procs[0];
            for( unsigned i = 1; i < ompNThreads; i++ ) {
                starting_node[i] = ending_node[i - 1];
                ending_node[i]   = starting_node[i] + npoint_procs[i];
            }

            int numflag, nparts, edgecut, wgtflag, ncon;

            int* vtxdist = new int[static_cast<unsigned>( ompNThreads ) + 1u];
            int* part    = new int[nPoint];

            real_t ubvec;
            real_t* tpwgts = new real_t[ompNThreads];

            wgtflag = 0;
            numflag = 0;
            ncon    = 1;
            ubvec   = static_cast<real_t>( 1.05 );
            nparts  = ompNThreads;
            int options[METIS_NOPTIONS];
            METIS_SetDefaultOptions( options );
            options[1] = 0;

            for( unsigned i = 0; i < ompNThreads; i++ ) {
                tpwgts[i] = static_cast<real_t>( 1.0 ) / static_cast<real_t>( ompNThreads );
            }

            vtxdist[0] = 0;
            for( unsigned i = 0; i < ompNThreads; i++ ) {
                vtxdist[i + 1] = static_cast<int>( ending_node[i] );
            }

            ParMETIS_V3_PartKway(
                vtxdist, xadj, adjncy, nullptr, nullptr, &wgtflag, &numflag, &ncon, &nparts, tpwgts, &ubvec, options, &edgecut, part, &comm );

            _colors.resize( nPoint );
            for( unsigned i = 0; i < nPoint; ++i ) {
                _colors[i] = static_cast<unsigned>( part[i] );
            }

            delete[] vtxdist;
            delete[] part;
            delete[] tpwgts;
        }
        else {    // use metis
            // METIS_PartGraphKway();
        }
    }
    else {
        _colors.resize( nPoint, 0u );
    }
}

unsigned Mesh::GetDim() const { return _dim; }
unsigned Mesh::GetNumCells() const { return _numCells; }
unsigned Mesh::GetNumNodes() const { return _numNodes; }
unsigned Mesh::GetNumNodesPerCell() const { return _numNodesPerCell; }

const std::vector<Vector>& Mesh::GetNodes() const { return _nodes; }
const std::vector<std::vector<unsigned>>& Mesh::GetCells() const { return _cells; }
const std::vector<double>& Mesh::GetCellAreas() const { return _cellAreas; };
