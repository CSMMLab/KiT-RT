#include "common/mesh.hpp"
#include "common/config.hpp"
#include "common/io.hpp"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <map>
#ifdef BUILD_MPI
#include <mpi.h>
#endif
#include <omp.h>
#include <set>

Mesh::Mesh( const Config* settings,
            std::vector<Vector> nodes,
            std::vector<std::vector<unsigned>> cells,
            std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> boundaries )
    : _dim( nodes[0].size() ), _numCells( cells.size() ), _numNodes( nodes.size() ), _numNodesPerCell( cells[0].size() ),
      _numBoundaries( boundaries.size() ), _ghostCellID( _numCells ), _nodes( nodes ), _cells( cells ), _boundaries( boundaries ) {

    _settings = settings;
    if( _dim == 2 ) {
        _numNodesPerBoundary = 2u;
    }
    else {
        ErrorMessages::Error( "Unsupported mesh dimension!", CURRENT_FUNCTION );
    }
    int nprocs = 1;
    int rank   = 0;
#ifdef BUILD_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif
    if( rank == 0 ) {

        auto log = spdlog::get( "event" );
        log->info( "| Compute cell areas..." );
    }
    ComputeCellAreas();
    if( rank == 0 ) {
        auto log = spdlog::get( "event" );
        log->info( "| Compute cell midpoints..." );
    }
    ComputeCellMidpoints();

    // Connectivity
    std::string connectivityFile = _settings->GetMeshFile();
    size_t lastDotIndex          = connectivityFile.find_last_of( '.' );
    connectivityFile             = connectivityFile.substr( 0, lastDotIndex );
    connectivityFile += ".con";
    if( !std::filesystem::exists( connectivityFile ) || _settings->GetForcedConnectivity() ) {
        if( rank == 0 ) {
            auto log = spdlog::get( "event" );

            log->info( "| Compute mesh connectivity..." );
        }
        ComputeConnectivity();    // Computes  _cellNeighbors, _cellInterfaceMidPoints, _cellNormals, _cellBoundaryTypes
        if( !_settings->GetForcedConnectivity() ) {
            if( rank == 0 ) {
                auto log = spdlog::get( "event" );

                log->info( "| Save mesh connectivity to file " + connectivityFile );
                WriteConnecitivityToFile(
                    connectivityFile, _cellNeighbors, _cellInterfaceMidPoints, _cellNormals, _cellBoundaryTypes, _numCells, _dim );
            }
        }
    }
    else {
        // Resize the outer vector to have nCells elements
        _cellNeighbors.resize( _numCells );
        _cellInterfaceMidPoints.resize( _numCells );
        _cellNormals.resize( _numCells );
        _cellBoundaryTypes.resize( _numCells );
        if( rank == 0 ) {
            auto log = spdlog::get( "event" );
            log->info( "| Load mesh connectivity from file " + connectivityFile );
        }
        LoadConnectivityFromFile(
            connectivityFile, _cellNeighbors, _cellInterfaceMidPoints, _cellNormals, _cellBoundaryTypes, _numCells, _numNodesPerCell, _dim );
    }
    if( rank == 0 ) {
        auto log = spdlog::get( "event" );
        log->info( "| Compute boundary..." );
    }
    ComputeBounds();
    if( rank == 0 ) {
        auto log = spdlog::get( "event" );
        log->info( "| Mesh created." );
    }
}

Mesh::~Mesh() {}

void Mesh::ComputeConnectivity() {
    int nprocs = 1;
    int rank   = 0;
#ifdef BUILD_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif

    unsigned comm_size = 1;    // No MPI implementation right now
    // determine number/chunk size and indices of cells treated by each mpi thread (deactivated for now)
    unsigned chunkSize    = _numCells;    // std::ceil( static_cast<float>( _numCells ) / static_cast<float>( comm_size ) );
    unsigned mpiCellStart = 0;            // comm_rank * chunkSize;
    unsigned mpiCellEnd   = _numCells;    // std::min( ( comm_rank + 1 ) * chunkSize, _numCells );

    // 'flat' vectors are a flattened representation of the neighbors numCells<numNodesPerCell> nested vectors; easier for MPI
    // 'part' vectors store information for each single MPI thread
    std::vector<int> neighborsFlatPart( _numNodesPerCell * chunkSize, -1 );
    std::vector<Vector> normalsFlatPart( _numNodesPerCell * chunkSize, Vector( _dim, -1.0 ) );
    std::vector<Vector> interfaceMidFlatPart( _numNodesPerCell * chunkSize, Vector( _dim, -1.0 ) );

    // ++++++
    std::vector<std::vector<int>> neighbors( _numCells );
    std::vector<std::vector<Vector>> normals( _numCells );
    std::vector<std::vector<Vector>> interfaceMid( _numCells );

    std::map<std::pair<unsigned, unsigned>, std::vector<unsigned>> edge_to_cells;
    //  Step 1: Build a mapping from edges to cells
    if( rank == 0 ) {
        auto log = spdlog::get( "event" );
        log->info( "| ...map edges to cells..." );
    }
    for( unsigned i = 0; i < _numCells; ++i ) {
        const auto& cell = _cells[i];
        for( unsigned j = 0; j < _numNodesPerCell; ++j ) {
            unsigned v1                        = cell[j];
            unsigned v2                        = cell[( j + 1 ) % _numNodesPerCell];
            std::pair<unsigned, unsigned> edge = std::minmax( v1, v2 );
            edge_to_cells[edge].push_back( i );
        }
    }

    // Step 2: Determine neighbors
    if( rank == 0 ) {
        auto log = spdlog::get( "event" );
        log->info( "| ...determine neighbors of cells..." );
    }

    for( const auto& item : edge_to_cells ) {
        const auto& cell_list     = item.second;
        const auto& nodes_of_edge = item.first;
        // std::cout << "nodes_of_edge: " << nodes_of_edge.first << " " << nodes_of_edge.second << std::endl;
        if( cell_list.size() == 2 ) {
            if( cell_list[0] == cell_list[1] ) {
                ErrorMessages::Error( "Error", CURRENT_FUNCTION );
            }
            unsigned cell1 = cell_list[0];
            unsigned cell2 = cell_list[1];

            neighbors[cell1].push_back( cell2 );
            neighbors[cell2].push_back( cell1 );

            normals[cell1].push_back(
                ComputeOutwardFacingNormal( _nodes[nodes_of_edge.first], _nodes[nodes_of_edge.second], _cellMidPoints[cell1] ) );
            normals[cell2].push_back(
                ComputeOutwardFacingNormal( _nodes[nodes_of_edge.first], _nodes[nodes_of_edge.second], _cellMidPoints[cell2] ) );
            interfaceMid[cell1].push_back( ComputeCellInterfaceMidpoints( _nodes[nodes_of_edge.first], _nodes[nodes_of_edge.second] ) );
            interfaceMid[cell2].push_back( ComputeCellInterfaceMidpoints( _nodes[nodes_of_edge.first], _nodes[nodes_of_edge.second] ) );
        }
        else if( /* condition */ cell_list.size() == 1 ) {

            unsigned cell1 = cell_list[0];
            neighbors[cell1].push_back( _ghostCellID );    // neighbor must be a ghost cell
            normals[cell1].push_back(
                ComputeOutwardFacingNormal( _nodes[nodes_of_edge.first], _nodes[nodes_of_edge.second], _cellMidPoints[cell1] ) );
            interfaceMid[cell1].push_back( ComputeCellInterfaceMidpoints( _nodes[nodes_of_edge.first], _nodes[nodes_of_edge.second] ) );
        }
        else {
            ErrorMessages::Error( "More than 2 cells share an edge: " + std::to_string( cell_list.size() ), CURRENT_FUNCTION );
        }
    }
    // check for any unassigned faces
    _cellNeighbors.resize( _numCells );
    for( unsigned i = 0; i < _numCells; ++i ) {
        _cellNeighbors[i].resize( _numNodesPerCell );
        if( neighbors[i].size() != _numNodesPerCell ) {
            ErrorMessages::Error( "Not " + std::to_string( _numNodesPerCell ) + " neighbors detected: " + std::to_string( i ), CURRENT_FUNCTION );
        }
        if( normals[i].size() != _numNodesPerCell ) {
            ErrorMessages::Error( "Not" + std::to_string( _numNodesPerCell ) + " normals detected: " + std::to_string( i ), CURRENT_FUNCTION );
        }
        if( interfaceMid[i].size() != _numNodesPerCell ) {
            ErrorMessages::Error( "Not" + std::to_string( _numNodesPerCell ) + " interfaceMid detected: " + std::to_string( i ), CURRENT_FUNCTION );
        }
        for( unsigned j = 0; j < _numNodesPerCell; ++j ) {
            if( neighbors[i][j] == -1 ) {
                ErrorMessages::Error( "Face not assigned: " + std::to_string( i ) + " " + std::to_string( j ), CURRENT_FUNCTION );
            }
            _cellNeighbors[i][j] = neighbors[i][j];
        }
    }
    _cellNormals            = normals;
    _cellInterfaceMidPoints = interfaceMid;

    /*
    // pre sort cells and boundaries; sorting is needed for std::set_intersection
    log->info( "| ...sort cells..." );

    auto sortedCells( _cells );
#pragma omp parallel for
    for( unsigned i = 0; i < _numCells; ++i ) {
        std::sort( sortedCells[i].begin(), sortedCells[i].end() );
    }
    std::vector<std::vector<unsigned>> sortedBoundaries;
    for( unsigned i = 0; i < _numBoundaries; ++i ) {
        sortedBoundaries.push_back( _boundaries[i].second );
        std::sort( sortedBoundaries[i].begin(), sortedBoundaries[i].end() );
    }

    // save which cell has which nodes
    log->info( "| ...connect cells to nodes..." );
    blaze::CompressedMatrix<bool> connMat( _numCells, _numNodes );

    // #pragma omp parallel for
    for( unsigned i = mpiCellStart; i < mpiCellEnd; ++i ) {
        for( auto j : _cells[i] ) connMat.set( i, j, true );
    }

    // determine neighbor cells and normals with MPI and OpenMP
    log->info( "| ...determine neighbors of cells..." );

#pragma omp parallel for schedule( guided )
    for( unsigned i = mpiCellStart; i < mpiCellEnd; ++i ) {
        std::vector<unsigned>* cellsI = &sortedCells[i];
        unsigned ctr                  = 0;
        for( unsigned j = 0; j < _numCells; ++j ) {
            if( i == j )
                continue;
            else if( ctr == _numNodesPerCell )
                break;
            else if( static_cast<unsigned>( blaze::dot( blaze::row( connMat, i ), blaze::row( connMat, j ) ) ) ==
                     _numNodesPerBoundary ) {    // in 2D cells are neighbors if they share two nodes std::vector<unsigned>* cellsJ =
                                                 // &sortedCells[j];
                // which are the two common nodes and which edge belongs to them?
                std::vector<unsigned>* cellsJ = &sortedCells[j];
                std::vector<unsigned> commonElements;    // vector of nodes that are shared by cells i and j
                commonElements.reserve( _numNodesPerBoundary );
                std::set_intersection( cellsI->begin(),
                                       cellsI->end(),
                                       cellsJ->begin(),
                                       cellsJ->end(),
                                       std::back_inserter( commonElements ) );    // find common nodes of two cells
                // determine unused index
                unsigned pos0 = _numNodesPerCell * ( i - mpiCellStart );
                unsigned pos  = pos0;
                while( neighborsFlatPart[pos] != -1 && pos < pos0 + _numNodesPerCell - 1 && pos < chunkSize * _numNodesPerCell - 1 ) {
                    pos++;    // neighbors should be at same edge position for cells i AND j
                }
                neighborsFlatPart[pos] = j;
                // compute normal vector
                normalsFlatPart[pos]      = ComputeOutwardFacingNormal( _nodes[commonElements[0]], _nodes[commonElements[1]], _cellMidPoints[i] );
                interfaceMidFlatPart[pos] = ComputeCellInterfaceMidpoints( _nodes[commonElements[0]], _nodes[commonElements[1]] );
                ctr++;
            }
        }

        // boundaries are treated similarly to normal cells, but need a special treatment due to the absence of a neighboring cell
        for( unsigned k = 0; k < _boundaries.size(); ++k ) {
            std::vector<unsigned>* bNodes = &sortedBoundaries[k];
            for( unsigned j = 0; j < _boundaries[k].second.size(); ++j ) {
                std::vector<unsigned> commonElements;    // vector of nodes that are shared by cells i and j
                commonElements.reserve( _numNodesPerBoundary );
                std::set_intersection( cellsI->begin(),
                                       cellsI->end(),
                                       bNodes->begin(),
                                       bNodes->end(),
                                       std::back_inserter( commonElements ) );    // find common nodes of two cells
                // _boundaries[k].second has all boundary nodes of boundary k. Therefore if all cell nodes lie on the boundary, the number of
                // common nodes can be 3 for triangles, 4 for quadrangles etc
                if( commonElements.size() >= _numNodesPerBoundary && commonElements.size() <= _numNodesPerCell ) {
                    unsigned pos0 = _numNodesPerCell * ( i - mpiCellStart );
                    unsigned pos  = pos0;
                    while( neighborsFlatPart[pos] != -1 && pos < pos0 + _numNodesPerCell - 1 && pos < chunkSize * _numNodesPerCell - 1 ) {
                        pos++;
                    }
                    neighborsFlatPart[pos]    = _ghostCellID;
                    normalsFlatPart[pos]      = ComputeOutwardFacingNormal( _nodes[commonElements[0]], _nodes[commonElements[1]],
_cellMidPoints[i] ); interfaceMidFlatPart[pos] = ComputeCellInterfaceMidpoints( _nodes[commonElements[0]], _nodes[commonElements[1]] );
                }
            }
        }
    }

    // gather distributed data on all MPI threads
    std::vector<int> neighborsFlat( _numNodesPerCell * chunkSize * comm_size, -1 );
    std::vector<Vector> normalsFlat( _numNodesPerCell * chunkSize * comm_size, Vector( _dim, 0.0 ) );
    std::vector<Vector> interfaceMidFlat( _numNodesPerCell * chunkSize * comm_size, Vector( _dim, 0.0 ) );
    // if( comm_size == 1 ) {    // can be done directly if there is only one MPI thread
    neighborsFlat.assign( neighborsFlatPart.begin(), neighborsFlatPart.end() );
    normalsFlat.assign( normalsFlatPart.begin(), normalsFlatPart.end() );
    interfaceMidFlat.assign( interfaceMidFlatPart.begin(), interfaceMidFlatPart.end() );
    //}
    // else {
    //    MPI_Allgather( neighborsFlatPart.data(),
    //                   _numNodesPerCell * chunkSize,
    //                   MPI_INT,
    //                   neighborsFlat.data(),
    //                   _numNodesPerCell * chunkSize,
    //                   MPI_INT,
    //                   MPI_COMM_WORLD );
    //}

    // check for any unassigned faces
    if( std::any_of( neighborsFlat.begin(), neighborsFlat.end(), []( int i ) { return i == -1; } ) ) {    // if any entry in neighborsFlat is -1
        for( unsigned idx = 0; idx < neighborsFlat.size(); ++idx ) {
            ErrorMessages::Error( "Detected unassigned faces at index " + std::to_string( idx ) + " !", CURRENT_FUNCTION );
        }
    }

    // reorder neighbors and normals into nested structure
    _cellNeighbors.resize( _numCells );
    _cellNormals.resize( _numCells );
    _cellInterfaceMidPoints.resize( _numCells );

    for( unsigned i = 0; i < neighborsFlat.size(); ++i ) {
        unsigned IDi = static_cast<unsigned>( i / static_cast<double>( _numNodesPerCell ) );
        unsigned IDj = neighborsFlat[i];
        if( IDi == IDj ) continue;     // avoid self assignment
        if( IDj == _ghostCellID ) {    // cell is boundary cell
            // if( std::find( _cellNeighbors[IDi].begin(), _cellNeighbors[IDi].end(), _ghostCellID ) == _cellNeighbors[IDi].end() ) {
            _cellNeighbors[IDi].push_back( _ghostCellID );
            _cellNormals[IDi].push_back( normalsFlat[i] );
            _cellInterfaceMidPoints[IDi].push_back( interfaceMidFlat[i] );
            // temp++;
            // }
            // std::cout << temp << "\n";
        }
        else {    // normal cell neighbor
            // if( std::find( _cellNeighbors[IDi].begin(), _cellNeighbors[IDi].end(), IDj ) == _cellNeighbors[IDi].end() ) {
            _cellNeighbors[IDi].push_back( neighborsFlat[i] );
            _cellNormals[IDi].push_back( normalsFlat[i] );
            _cellInterfaceMidPoints[IDi].push_back( interfaceMidFlat[i] );
            //}
        }
    }
*/
    // assign boundary types to all cells
    _cellBoundaryTypes.resize( _numCells, BOUNDARY_TYPE::NONE );
#pragma omp parallel for
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
#pragma omp parallel for
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
#pragma omp parallel for
    for( unsigned j = 0; j < _numCells; ++j ) {               // loop over cells
        for( unsigned l = 0; l < _cells[j].size(); ++l ) {    // loop over nodes of the cell
            _cellMidPoints[j] = _cellMidPoints[j] + _nodes[_cells[j][l]];
        }
        _cellMidPoints[j] = _cellMidPoints[j] / static_cast<double>( _cells[j].size() );    // arithmetic mean of node coord
    }
}

Vector Mesh::ComputeCellInterfaceMidpoints( const Vector& nodeA, const Vector& nodeB ) {
    Vector interfaceMidPt( _dim, 0.0 );
    interfaceMidPt = 0.5 * ( nodeA + nodeB );
    return interfaceMidPt;
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
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _numCells; ++idx_cell ) {
        for( unsigned idx_sys = 0; idx_sys < nq; ++idx_sys ) {
            psiDerX[idx_cell][idx_sys] = 0.0;
            psiDerY[idx_cell][idx_sys] = 0.0;
            if( _cellBoundaryTypes[idx_cell] != 2 ) continue;    // skip ghost cells
            //  compute derivative by summing over cell boundary using Green Gauss Theorem
            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); ++idx_nbr ) {
                psiDerX[idx_cell][idx_sys] +=
                    0.5 * ( psi[idx_cell][idx_sys] + psi[_cellNeighbors[idx_cell][idx_nbr]][idx_sys] ) * _cellNormals[idx_cell][idx_nbr][0];
                psiDerY[idx_cell][idx_sys] +=
                    0.5 * ( psi[idx_cell][idx_sys] + psi[_cellNeighbors[idx_cell][idx_nbr]][idx_sys] ) * _cellNormals[idx_cell][idx_nbr][1];
            }
            psiDerX[idx_cell][idx_sys] /= _cellAreas[idx_cell];
            psiDerY[idx_cell][idx_sys] /= _cellAreas[idx_cell];
        }
    }
}

void Mesh::ComputeSlopes1D( unsigned nq, VectorVector& psiDerX, const VectorVector& psi ) const {
    // assume equidistant ordered mesh
    double dx = _cellMidPoints[2][0] - _cellMidPoints[1][0];

#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _numCells; ++idx_cell ) {
        for( unsigned idx_sys = 0; idx_sys < nq; ++idx_sys ) {
            psiDerX[idx_cell][idx_sys] = 0.0;
            if( _cellNeighbors[idx_cell].size() <= 2 ) {    // right neighbor + ghostcell ==> its a boundary cell
                continue;                                   // skip computation
            }
            // compute derivative by second order difference
            psiDerX[idx_cell][idx_sys] = ( psi[idx_cell + 1][idx_sys] - 2 * psi[idx_cell][idx_sys] + psi[idx_cell - 1][idx_sys] ) / ( dx * dx );
        }
    }
}

void Mesh::ComputeLimiter(
    unsigned nSys, const VectorVector& solDx, const VectorVector& solDy, const VectorVector& sol, VectorVector& limiter ) const {
    double const eps = 1e-10;
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _numCells; idx_cell++ ) {
        if( _cellBoundaryTypes[idx_cell] != 2 ) {
            for( unsigned idx_sys = 0; idx_sys < nSys; idx_sys++ ) {
                limiter[idx_cell][idx_sys] = 0.0;    // turn to first order on boundaries
            }
            continue;    // skip computation
        }

        for( unsigned idx_sys = 0; idx_sys < nSys; idx_sys++ ) {
            double r      = 0.0;
            double minSol = sol[idx_cell][idx_sys];
            double maxSol = sol[idx_cell][idx_sys];

            Vector localLimiter( _numNodesPerCell, 0.0 );
            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); idx_nbr++ ) {
                // Compute ptswise max and minimum solultion values of current and neighbor cells
                unsigned glob_nbr = _cellNeighbors[idx_cell][idx_nbr];
                if( sol[glob_nbr][idx_sys] > maxSol ) {
                    maxSol = sol[glob_nbr][idx_sys];
                }
                if( sol[glob_nbr][idx_sys] < minSol ) {
                    minSol = sol[glob_nbr][idx_sys];
                }
            }

            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); idx_nbr++ ) {
                // Compute value at interface midpoint, called gaussPt
                double gaussPt = 0.0;
                // gauss point is at cell vertex
                gaussPt = ( solDx[idx_cell][idx_sys] * ( _nodes[_cells[idx_cell][idx_nbr]][0] - _cellMidPoints[idx_cell][0] ) +
                            solDy[idx_cell][idx_sys] * ( _nodes[_cells[idx_cell][idx_nbr]][1] - _cellMidPoints[idx_cell][1] ) );

                // Compute limiter input
                if( std::abs( gaussPt ) > eps ) {
                    if( gaussPt > 0.0 ) {
                        r = ( maxSol - sol[idx_cell][idx_sys] ) / gaussPt;
                    }
                    else {
                        r = ( minSol - sol[idx_cell][idx_sys] ) / gaussPt;
                    }
                }
                else {
                    r = 1.0;
                }
                if( r < 0.0 ) {
                    std::cout << "r <0.0 \n";    // if this happens there is a bug or a deformend mesh
                }
                localLimiter[idx_nbr] = std::min( r, 1.0 );    // LimiterBarthJespersen( r );
            }
            // get smallest limiter
            limiter[idx_cell][idx_sys] = localLimiter[0];
            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); idx_nbr++ ) {
                if( localLimiter[idx_nbr] < limiter[idx_cell][idx_sys] ) limiter[idx_cell][idx_sys] = localLimiter[idx_nbr];
            }
        }
    }
}

void Mesh::ComputeLimiter1D( unsigned nSys, const VectorVector& sol, VectorVector& limiter ) const {
    double const eps = 1e-10;
#pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _numCells; idx_cell++ ) {
        for( unsigned idx_sys = 0; idx_sys < nSys; idx_sys++ ) {
            double r = 0.0;
            if( _cellNeighbors[idx_cell].size() <= 2 ) {    // right neighbor + ghostcell ==> its a boundary cell
                limiter[idx_cell][idx_sys] = 0.0;           // turn to first order on boundaries
                continue;                                   // skip computation
            }
            double up   = sol[idx_cell][idx_sys] - sol[_cellNeighbors[idx_cell][0]][idx_sys];
            double down = sol[_cellNeighbors[idx_cell][1]][idx_sys] - sol[idx_cell][idx_sys];

            up > 0 ? up += eps : up -= eps;          // to prevent divbyzero
            down > 0 ? down += eps : down -= eps;    // to prevent divbyzero

            r                          = up / down;
            limiter[idx_cell][idx_sys] = std::max( std::min( r, 1.0 ), 0.0 );    // minmod limiter
        }
    }
}

// double LimiterBarthJespersen( double r ) { return std::min( r, 1.0 ); }

void Mesh::ComputeBounds() {
    _bounds = std::vector( _dim, std::make_pair( 0.0, 1.0 ) );    // Currently not used
    //_bounds = std::vector( _dim, std::make_pair( std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() ) );
    //
    //     for( unsigned i = 0; i < _numNodes; ++i ) {
    //         for( unsigned j = 0; j < _dim; ++j ) {
    //             if( _nodes[i][j] < _bounds[j].first ) _bounds[j].first = _nodes[i][j];
    //             if( _nodes[i][j] > _bounds[j].second ) _bounds[j].second = _nodes[i][j];
    //         }
    //     }
}

const std::vector<Vector>& Mesh::GetNodes() const { return _nodes; }
const std::vector<Vector>& Mesh::GetCellMidPoints() const { return _cellMidPoints; }
const std::vector<std::vector<unsigned>>& Mesh::GetCells() const { return _cells; }
const std::vector<double>& Mesh::GetCellAreas() const { return _cellAreas; }
const std::vector<std::vector<unsigned>>& Mesh::GetNeighbours() const { return _cellNeighbors; }
const std::vector<std::vector<Vector>>& Mesh::GetNormals() const { return _cellNormals; }
const std::vector<BOUNDARY_TYPE>& Mesh::GetBoundaryTypes() const { return _cellBoundaryTypes; }
const std::vector<std::pair<double, double>> Mesh::GetBounds() const { return _bounds; }
const std::vector<std::vector<Vector>> Mesh::GetInterfaceMidPoints() const { return _cellInterfaceMidPoints; }

void Mesh::SetBoundaryType( int idx_cell, BOUNDARY_TYPE boundary_type ) { _cellBoundaryTypes[idx_cell] = boundary_type; }

double Mesh::GetDistanceToOrigin( unsigned idx_cell ) const {
    double distance = 0.0;
    for( unsigned idx_dim = 0; idx_dim < _dim; idx_dim++ ) {
        distance += _cellMidPoints[idx_cell][idx_dim] * _cellMidPoints[idx_cell][idx_dim];
    }
    return sqrt( distance );
}


unsigned Mesh::GetCellOfKoordinate( const double x, const double y ) const {
    // Experimental parallel implementation
    unsigned koordinate_cell_id = std::numeric_limits<unsigned>::max();
    bool found                  = false;

//#pragma omp parallel for shared( found )
    for( unsigned idx_cell = 0; idx_cell < _numCells; idx_cell++ ) {
        if( IsPointInsideCell( idx_cell, x, y ) ) {
            //#pragma omp critical
            {
                if( !found ) {
                    koordinate_cell_id = idx_cell;
                    found              = true;
                }
            }
        }
        // Check if cancellation has been requested
    }
    if( !found ) {
        ErrorMessages::Error( "Probing point (" + std::to_string( x ) + "," + std::to_string( y ) + ") is not contained in mesh.", CURRENT_FUNCTION );
    }
    return koordinate_cell_id;
}

std::vector<unsigned> Mesh::GetCellsofBall( const double x, const double y, const double r ) const {
    std::vector<unsigned> cells_in_ball;

    // Experimental parallel implementation
    // #pragma omp parallel for
    for( unsigned idx_cell = 0; idx_cell < _numCells; idx_cell++ ) {
        // Assume GetCellCenter returns the center coordinates of the cell
        double cell_x = _cellMidPoints[idx_cell][0];
        double cell_y = _cellMidPoints[idx_cell][1];

        // Calculate the distance from the cell center to the point (x, y)
        double distance = std::sqrt( ( cell_x - x ) * ( cell_x - x ) + ( cell_y - y ) * ( cell_y - y ) );

        if( distance <= r ) {
            // #pragma omp critical
            { cells_in_ball.push_back( idx_cell ); }
        }
    }

    if( cells_in_ball.empty() ) {    // take the only cell that contains the point
        std::cout << "No cells found within the ball centered at (" << x << "," << y << ") with radius " << r << "." << std::endl;
        std::cout << "Taking the only cell that contains the point." << std::endl;
        cells_in_ball.push_back( GetCellOfKoordinate( x, y ) );
        // ErrorMessages::Error( "No cells found within the ball centered at (" + std::to_string( x ) + "," + std::to_string( y ) + ") with radius " +
        // std::to_string( r ) + ".",
        // CURRENT_FUNCTION );
    }

    return cells_in_ball;
}

bool Mesh::IsPointInsideCell( unsigned idx_cell, double x, double y ) const {
    bool inside = false;
    if( _numNodesPerCell == 3 ) {
        inside = isPointInTriangle( x,
                                    y,
                                    _nodes[_cells[idx_cell][0]][0],
                                    _nodes[_cells[idx_cell][0]][1],
                                    _nodes[_cells[idx_cell][1]][0],
                                    _nodes[_cells[idx_cell][1]][1],
                                    _nodes[_cells[idx_cell][2]][0],
                                    _nodes[_cells[idx_cell][2]][1] );
    }
    else if( _numNodesPerCell == 4 ) {    // quadrilateral cell  +> divide into 2 triangles
        inside = isPointInTriangle( x,
                                    y,
                                    _nodes[_cells[idx_cell][0]][0],
                                    _nodes[_cells[idx_cell][0]][1],
                                    _nodes[_cells[idx_cell][1]][0],
                                    _nodes[_cells[idx_cell][1]][1],
                                    _nodes[_cells[idx_cell][2]][0],
                                    _nodes[_cells[idx_cell][2]][1] ) ||
                 isPointInTriangle( x,
                                    y,
                                    _nodes[_cells[idx_cell][0]][0],
                                    _nodes[_cells[idx_cell][0]][1],
                                    _nodes[_cells[idx_cell][2]][0],
                                    _nodes[_cells[idx_cell][2]][1],
                                    _nodes[_cells[idx_cell][3]][0],
                                    _nodes[_cells[idx_cell][3]][1] );
    }
    else {
        ErrorMessages::Error( "Unsupported number of nodes per cell: " + std::to_string( _numNodesPerCell ), CURRENT_FUNCTION );
    }
    // if( inside ) {
    //     std::cout << "Cells: " << idx_cell << "pt: " << _cellMidPoints[idx_cell][0] << " " << _cellMidPoints[idx_cell][1] << " pt: " << x << " " <<
    //     y
    //               << "\n";
    // }

    return inside;
}

bool Mesh::isPointInTriangle( double x, double y, double x1, double y1, double x2, double y2, double x3, double y3 ) const {
    // Calculate the area of the triangle
    double d1, d2, d3;
    bool hasNeg, hasPos;

    d1 = ( x - x2 ) * ( y1 - y2 ) - ( x1 - x2 ) * ( y - y2 );
    d2 = ( x - x3 ) * ( y2 - y3 ) - ( x2 - x3 ) * ( y - y3 );
    d3 = ( x - x1 ) * ( y3 - y1 ) - ( x3 - x1 ) * ( y - y1 );

    hasNeg = ( d1 < 0 ) || ( d2 < 0 ) || ( d3 < 0 );
    hasPos = ( d1 > 0 ) || ( d2 > 0 ) || ( d3 > 0 );

    return !( hasNeg && hasPos );
}
