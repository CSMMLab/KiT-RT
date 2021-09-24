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
    // ComputeCellInterfaceMidpoints();
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
    std::vector<Vector> interfaceMidFlatPart( _numNodesPerCell * chunkSize, Vector( _dim, -1.0 ) );

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

    // save which cell has which nodes
    blaze::CompressedMatrix<bool> connMat( _numCells, _numNodes );
    for( unsigned i = mpiCellStart; i < mpiCellEnd; ++i ) {
        for( auto j : _cells[i] ) connMat.set( i, j, true );
    }

// determine neighbor cells and normals with MPI and OpenMP
#pragma omp parallel for
    for( unsigned i = mpiCellStart; i < mpiCellEnd; ++i ) {
        std::vector<unsigned>* cellsI = &sortedCells[i];
        unsigned ctr                  = 0;
        for( unsigned j = 0; j < _numCells; ++j ) {
            if( i == j )
                continue;
            else if( ctr == _numNodesPerCell )
                break;
            else if( static_cast<unsigned>( blaze::dot( blaze::row( connMat, i ), blaze::row( connMat, j ) ) ) ==
                     _numNodesPerBoundary ) {    // in 2D cells are neighbors if they share two nodes std::vector<unsigned>* cellsJ = &sortedCells[j];
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
                while( neighborsFlatPart[pos] != -1 && pos < pos0 + _numNodesPerCell - 1 && pos < chunkSize * _numNodesPerCell - 1 )
                    pos++;    // neighbors should be at same edge position for cells i AND j
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
                // _boundaries[k].second has all boundary nodes of boundary k. Therefore if all cell nodes lie on the boundary, the number of common
                // nodes can be 3 for triangles, 4 for quadrangles etc
                if( commonElements.size() >= _numNodesPerBoundary && commonElements.size() <= _numNodesPerCell ) {
                    unsigned pos0 = _numNodesPerCell * ( i - mpiCellStart );
                    unsigned pos  = pos0;
                    while( neighborsFlatPart[pos] != -1 && pos < pos0 + _numNodesPerCell - 1 && pos < chunkSize * _numNodesPerCell - 1 ) pos++;
                    neighborsFlatPart[pos]    = _ghostCellID;
                    normalsFlatPart[pos]      = ComputeOutwardFacingNormal( _nodes[commonElements[0]], _nodes[commonElements[1]], _cellMidPoints[i] );
                    interfaceMidFlatPart[pos] = ComputeCellInterfaceMidpoints( _nodes[commonElements[0]], _nodes[commonElements[1]] );
                }
            }
        }
    }

    // gather distributed data on all MPI threads
    std::vector<int> neighborsFlat( _numNodesPerCell * chunkSize * comm_size, -1 );
    std::vector<Vector> normalsFlat( _numNodesPerCell * chunkSize * comm_size, Vector( _dim, 0.0 ) );
    std::vector<Vector> interfaceMidFlat( _numNodesPerCell * chunkSize * comm_size, Vector( _dim, 0.0 ) );
    if( comm_size == 1 ) {    // can be done directly if there is only one MPI thread
        neighborsFlat.assign( neighborsFlatPart.begin(), neighborsFlatPart.end() );
        normalsFlat.assign( normalsFlatPart.begin(), normalsFlatPart.end() );
        interfaceMidFlat.assign( interfaceMidFlatPart.begin(), interfaceMidFlatPart.end() );
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
            if( std::find( _cellNeighbors[IDi].begin(), _cellNeighbors[IDi].end(), _ghostCellID ) == _cellNeighbors[IDi].end() ) {
                _cellNeighbors[IDi].push_back( _ghostCellID );
                _cellNormals[IDi].push_back( normalsFlat[i] );
                _cellInterfaceMidPoints[IDi].push_back( interfaceMidFlat[i] );
            }
        }
        else {    // normal cell neighbor
            if( std::find( _cellNeighbors[IDi].begin(), _cellNeighbors[IDi].end(), IDj ) == _cellNeighbors[IDi].end() ) {
                _cellNeighbors[IDi].push_back( IDj );
                _cellNormals[IDi].push_back( normalsFlat[i] );
                _cellInterfaceMidPoints[IDi].push_back( interfaceMidFlat[i] );
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
    for( unsigned idx_sys = 0; idx_sys < nq; ++idx_sys ) {
        for( unsigned idx_cell = 0; idx_cell < _numCells; ++idx_cell ) {
            psiDerX[idx_cell][idx_sys] = 0.0;
            psiDerY[idx_cell][idx_sys] = 0.0;

            // if( cell->IsBoundaryCell() ) continue; // skip ghost cells
            if( _cellBoundaryTypes[idx_cell] != 2 ) continue;    // skip ghost cells
            // compute derivative by summing over cell boundary
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

void Mesh::ComputeLimiter(
    unsigned nSys, const VectorVector& solDx, const VectorVector& solDy, const VectorVector& sol, VectorVector& limiter ) const {
    double r    = 0.0;
    double eps  = 1e-7;
    double sign = 0.0;
    for( unsigned idx_cell = 0; idx_cell < _numCells; idx_cell++ ) {
        for( unsigned idx_sys = 0; idx_sys < nSys; idx_sys++ ) {
            if( _cellBoundaryTypes[idx_cell] != 2 ) {
                limiter[idx_cell][idx_sys] = 0.0;    // turn to first order on boundaries
                continue;                            // skip computation
            }
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
                if( idx_cell == 1049 && idx_nbr == 1 ) {
                    // std::cout << "here\n";
                }
                double gaussPt = 0.0;
                // double dy      = solDy[idx_cell][idx_sys];
                // double dx      = solDx[idx_cell][idx_sys];
                // double rijx    = _cellInterfaceMidPoints[idx_cell][idx_nbr][0];
                // double rijy    = _cellInterfaceMidPoints[idx_cell][idx_nbr][1];
                // double cmx     = _cellMidPoints[idx_cell][0];
                // double cmy     = _cellMidPoints[idx_cell][1];
                // double curSol  = sol[idx_cell][idx_sys];
                //_nodes[_cells[idx_cell][idx_nbr]]
                // _interfaceMidPoints[idx_cell][idx_nbr]
                // gauss point is at cell vertex
                gaussPt = solDx[idx_cell][idx_sys] * ( _nodes[_cells[idx_cell][idx_nbr]][0] - _cellMidPoints[idx_cell][0] ) +
                          solDy[idx_cell][idx_sys] * ( _nodes[_cells[idx_cell][idx_nbr]][1] - _cellMidPoints[idx_cell][1] );
                // Compute limiter input
                if( gaussPt < 0.0 )
                    sign = -1.0;
                else
                    sign = 1.0;

                if( std::abs( gaussPt ) > eps ) {
                    if( gaussPt > 0.0 ) {
                        r = ( maxSol - sol[idx_cell][idx_sys] ) / gaussPt;
                    }
                    else if( gaussPt < 0.0 ) {
                        r = ( minSol - sol[idx_cell][idx_sys] ) / gaussPt;
                    }
                }
                else {
                    r = 1.0;
                }
                // r = 1.0;
                if( r < 0.0 ) {
                    std::cout << "r <0.0 \n";
                }
                localLimiter[idx_nbr] = std::min( r, 1.0 );    // LimiterBarthJespersen( r );
                // double epsVenka = ( 1 * sqrt( _cellAreas[idx_cell] ) );
                // epsVenka        = epsVenka * epsVenka * epsVenka;
                // double dMax     = maxSol - sol[idx_cell][idx_sys];
                // double dMin     = minSol - sol[idx_cell][idx_sys];
                //// venkat limiter
                // if( gaussPt > 0.0 ) {
                //    localLimiter[idx_nbr] = ( 1 / gaussPt ) * ( ( dMax * dMax + epsVenka * epsVenka ) * gaussPt + 2 * gaussPt * gaussPt * dMax ) /
                //                            ( dMax * dMax + 2 * gaussPt * gaussPt + dMax * gaussPt + epsVenka * epsVenka );
                //}
                // else if( gaussPt < 0.0 ) {
                //    localLimiter[idx_nbr] = ( 1 / gaussPt ) * ( ( dMin * dMin + epsVenka * epsVenka ) * gaussPt + 2 * gaussPt * gaussPt * dMin ) /
                //                            ( dMin * dMin + 2 * gaussPt * gaussPt + dMin * gaussPt + epsVenka * epsVenka );
                //}
                // else {
                //    localLimiter[idx_nbr] = 1.0;
                //}
            }
            // get smallest limiter
            limiter[idx_cell][idx_sys] = localLimiter[0];
            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); idx_nbr++ ) {
                if( localLimiter[idx_nbr] < limiter[idx_cell][idx_sys] ) limiter[idx_cell][idx_sys] = localLimiter[idx_nbr];
            }
            // check maximum principle
            for( unsigned idx_nbr = 0; idx_nbr < _cellNeighbors[idx_cell].size(); idx_nbr++ ) {
                double currLim = limiter[idx_cell][idx_sys];
                double dy      = solDy[idx_cell][idx_sys];
                double dx      = solDx[idx_cell][idx_sys];
                double rijx    = _cellInterfaceMidPoints[idx_cell][idx_nbr][0];
                double rijy    = _cellInterfaceMidPoints[idx_cell][idx_nbr][1];
                double cmx     = _cellMidPoints[idx_cell][0];
                double cmy     = _cellMidPoints[idx_cell][1];
                double curSol  = sol[idx_cell][idx_sys];
                double gaussPt = solDx[idx_cell][idx_sys] * ( _cellInterfaceMidPoints[idx_cell][idx_nbr][0] - _cellMidPoints[idx_cell][0] ) +
                                 solDy[idx_cell][idx_sys] * ( _cellInterfaceMidPoints[idx_cell][idx_nbr][1] - _cellMidPoints[idx_cell][1] );

                double psiL  = sol[idx_cell][idx_sys] + currLim * gaussPt;
                double psiL2 = curSol + currLim * ( dx * ( rijx - cmx ) + dy * ( rijy - cmy ) );

                if( psiL > maxSol ) {
                    // std::cout << "max principle hurt\n";
                    // gaussPt = solDx[idx_cell][idx_sys] * ( _nodes[_cells[idx_cell][idx_nbr]][0] - _cellMidPoints[idx_cell][0] ) +
                    //           solDy[idx_cell][idx_sys] * ( _nodes[_cells[idx_cell][idx_nbr]][1] - _cellMidPoints[idx_cell][1] );
                    // std::cout << "gaussPt" << gaussPt << "\n";
                    // std::cout << "enumMax" << maxSol - sol[idx_cell][idx_sys] << "\n";
                    // std::cout << "enumMin" << minSol - sol[idx_cell][idx_sys] << "\n";
                    // std::cout << "minSol" << minSol << "psiL" << psiL << "maxSol" << maxSol << "\n";
                    // limiter[idx_cell][idx_sys] = 0.0;
                }
                if( psiL < minSol ) {
                    // std::cout << "min principle hurt\n";
                    // std::cout << "gaussPt" << gaussPt << "\n";
                    // std::cout << "enumMax" << maxSol - sol[idx_cell][idx_sys] << "\n";
                    // std::cout << "enumMin" << minSol - sol[idx_cell][idx_sys] << "\n";
                    // std::cout << "minSol" << minSol << "psiL" << psiL << "maxSol" << maxSol << "\n";
                    // limiter[idx_cell][idx_sys] = 0.0;
                }
            }
        }
    }
}

// double LimiterBarthJespersen( double r ) { return std::min( r, 1.0 ); }

void Mesh::LimitSlopes( unsigned nSys, VectorVector& solDx, VectorVector& solDy, const VectorVector& sol ) const {
    /* Requires correct derivative values */
    double phi;
    double eps = 1e-10;    // safety epsilon, hard coded

    Vector deltaSolMax( nSys, 0.0 );
    Vector deltaSolMin( nSys, 0.0 );
    VectorVector gaussPt( nSys, Vector( _numNodesPerCell, 0.0 ) );
    VectorVector limiter( nSys, Vector( _numNodesPerCell, 0.0 ) );

    for( unsigned i = 0; i < _numCells; ++i ) {
        // skip boundary cells
        if( _cellBoundaryTypes[i] != 2 ) continue;
        // reset local values
        for( unsigned q = 0; q < nSys; ++q ) {
            deltaSolMax[q] = 0.0;
            deltaSolMin[q] = 0.0;
            for( unsigned j = 0; j < _numNodesPerCell; ++j ) {
                gaussPt[q][j] = 0.0;
                limiter[q][j] = 0.0;
            }
        }
        // calculate  largest difference of values between current and neighboring cells
        for( unsigned q = 0; q < nSys; ++q ) {
            for( unsigned j = 0; j < _numNodesPerCell; ++j ) {
                double deltaNbr = sol[_cellNeighbors[i][j]][q] - sol[i][q];
                if( deltaNbr > deltaSolMax[q] ) deltaSolMax[q] = deltaNbr;    // largest positive difference
                if( deltaNbr < deltaSolMin[q] ) deltaSolMin[q] = deltaNbr;    // largest negative difference
            }
        }
        // compute unconstrained value at each gauss point
        for( unsigned q = 0; q < nSys; ++q ) {
            for( unsigned j = 0; j < _numNodesPerCell; ++j ) {

                // step 2: 1st order gauss point
                gaussPt[q][j] = 0.5 * solDx[i][q] * ( _cellMidPoints[_cellNeighbors[i][j]][0] - _cellMidPoints[i][0] ) +
                                0.5 * solDy[i][q] * ( _cellMidPoints[_cellNeighbors[i][j]][1] - _cellMidPoints[i][1] );
            }
        }
        // compute venkatakrishnan limiter function
        for( unsigned q = 0; q < nSys; ++q ) {
            phi = std::numeric_limits<double>::infinity();
            for( unsigned j = 0; j < _numNodesPerCell; ++j ) {
                double y = 0.0;
                if( gaussPt[q][j] > 0.0 ) y = deltaSolMax[q] / gaussPt[q][j];
                if( gaussPt[q][j] < 0.0 )
                    y = deltaSolMin[q] / gaussPt[q][j];
                else
                    y = 1.0;
                limiter[q][j] = ( y * y + 2 * y ) / ( y * y + y + 2 );    // original venkatakrishnan
                // get minimum limiter
                if( limiter[q][j] < phi ) phi = limiter[q][j];
            }
            // step 5: limit the slope reconstructed from Gauss theorem
            solDx[i][q] *= phi;
            solDy[i][q] *= phi;

            for( unsigned j = 0; j < _numNodesPerCell; ++j ) {
                // Check if reconstruction is TVD
                // unsigned nbr_glob = _cellNeighbors[i][j];

                // double psiL = sol[i][q] + solDx[i][q] * ( _interfaceMidPoints[i][j][0] - _cellMidPoints[i][0] ) +
                //              solDy[i][q] * ( _interfaceMidPoints[i][j][1] - _cellMidPoints[i][1] );
                // double psiR = sol[nbr_glob][q] + solDx[nbr_glob][q] * ( _interfaceMidPoints[i][j][0] - _cellMidPoints[nbr_glob][0] ) +
                //              solDy[nbr_glob][q] * ( _interfaceMidPoints[i][j][1] - _cellMidPoints[nbr_glob][1] );

                // get max differnce between nbrs and recons

                double deltaSolMax_test  = sol[i][q];
                double deltaSolMin_test  = sol[i][q];
                double deltaSolMax_testR = sol[i][q];
                double deltaSolMin_testR = sol[i][q];

                for( unsigned l = 0; l < _numNodesPerCell; ++l ) {
                    double deltaNbrL = sol[_cellNeighbors[i][l]][q];
                    double deltaNbrR = sol[_cellNeighbors[i][l]][q];
                    if( deltaNbrL > deltaSolMax_test ) deltaSolMax_test = deltaNbrL;      // largest positive difference
                    if( deltaNbrL < deltaSolMin_test ) deltaSolMin_test = deltaNbrL;      // largest negative difference
                    if( deltaNbrR > deltaSolMax_testR ) deltaSolMax_testR = deltaNbrR;    // largest positive difference
                    if( deltaNbrR < deltaSolMin_testR ) deltaSolMin_testR = deltaNbrR;    // largest negative difference
                }
                // compare extrema
                // if( psiL > deltaSolMax_test + eps ) {
                //    std::cout << "tvd not ok psiL max\n";
                //    deltaSolMax_test = psiL;
                //}    // largest positive difference
                // if( psiL + eps < deltaSolMin_test ) {
                //    std::cout << "tvd not ok psiL min\n";
                //    deltaSolMin_test = psiL;
                //}    // largest negative difference
                // if( psiR > deltaSolMax_testR + eps ) {
                //    std::cout << "tvd not ok psiR max\n";
                //    deltaSolMax_testR = psiR;
                //}    // largest positive difference
                // if( psiR + eps < deltaSolMin_testR ) {
                //    std::cout << "tvd not ok psiR min\n";
                //    deltaSolMin_testR = psiR;
                //}    // largest negative difference
            }
        }
    }
}

/*
void Mesh::ReconstructSlopesU2( unsigned nSys, VectorVector& psiDerX, VectorVector& psiDerY, const VectorVector& psi ) const {

    double phi;
    double eps = 1e-3;

    // VectorVector dPsiMax = std::vector( _numCells, Vector( nq, 0.0 ) );
    // VectorVector dPsiMin = std::vector( _numCells, Vector( nq, 0.0 ) );
    // std::vector<std::vector<Vector>> psiSample( _numCells, std::vector<Vector>( nSys, Vector( _numNodesPerCell, 0.0 ) ) );
    // std::vector<std::vector<Vector>> phiSample( _numCells, std::vector<Vector>( nSys, Vector( _numNodesPerCell, 0.0 ) ) );

    for( unsigned i = 0; i < _numCells; ++i ) {
        // set local values
        Vector dPsiMax( nSys, 0.0 );
        Vector dPsiMin( nSys, 0.0 );
        VectorVector delta( nSys, Vector( _numNodesPerCell, 0.0 ) );
        VectorVector limiter( nSys, Vector( _numNodesPerCell, 0.0 ) );
        // reset derivatives
        for( unsigned q = 0; q < nSys; ++q ) {
            psiDerX[i][q] = 0.0;
            psiDerY[i][q] = 0.0;
        }
        // skip boundary cells
        if( _cellBoundaryTypes[i] != 2 ) continue;
        // calculate limiter
        for( unsigned q = 0; q < nSys; ++q ) {
            // version 2: Venkatakrishnan limiter
            // step 1: calculate psi difference around neighbors and theoretical derivatives by Gauss theorem
            for( unsigned j = 0; j < _cellNeighbors[i].size(); ++j ) {
                if( psi[_cellNeighbors[i][j]][q] > dPsiMax[q] ) dPsiMax[q] = psi[_cellNeighbors[i][j]][q];
                if( psi[_cellNeighbors[i][j]][q] < dPsiMin[q] ) dPsiMin[q] = psi[_cellNeighbors[i][j]][q];
                dPsiMax[q] -= -psi[i][q];
                dPsiMin[q] -= -psi[i][q];

                psiDerX[i][q] += 0.5 * ( psi[i][q] + psi[_cellNeighbors[i][j]][q] ) * _cellNormals[i][j][0] / _cellAreas[i];
                psiDerY[i][q] += 0.5 * ( psi[i][q] + psi[_cellNeighbors[i][j]][q] ) * _cellNormals[i][j][1] / _cellAreas[i];
            }

            for( unsigned j = 0; j < _cellNeighbors[i].size(); ++j ) {

                // step 2: choose sample points Delta
                delta[q][j] = psiDerX[i][q] * ( _interfaceMidPoints[i][j][0] - _cellMidPoints[i][0] ) +
                              psiDerY[i][q] * ( _interfaceMidPoints[i][j][1] - _cellMidPoints[i][1] );

                // step 3: calculate Phi_ij at sample points
                if( delta[q][j] > 0.0 ) {
                    limiter[q][j] = ( dPsiMax[q] * dPsiMax[q] + 2.0 * dPsiMax[q] * delta[q][j] + eps ) /
                                    ( dPsiMax[q] * dPsiMax[q] + 2.0 * delta[q][j] * delta[q][j] + dPsiMax[q] * delta[q][j] + eps );
                }
                else if( delta[q][j] < 0.0 ) {
                    limiter[q][j] = ( dPsiMin[q] * dPsiMin[q] + 2.0 * dPsiMin[q] * delta[q][j] + eps ) /
                                    ( dPsiMin[q] * dPsiMin[q] + 2.0 * delta[q][j] * delta[q][j] + dPsiMin[q] * delta[q][j] + eps );
                }
                else {
                    limiter[q][j] = 1.0;
                }
                std::cout << limiter[q][j] << std::endl;

                if( limiter[q][j] < 0 ) {
                    std::cout << "error:" << limiter[q][j] << std::endl;
                }
                if( limiter[q][j] < -1.0 ) {
                    std::cout << "error:" << limiter[q][j] << std::endl;
                }
            }

            // step 4: find minimum limiter function phi
            phi = min( limiter[q] );
            // phi = fmin( 0.5, abs(phi) );

            // step 5: limit the slope reconstructed from Gauss theorem
            psiDerX[i][q] *= phi;
            psiDerY[i][q] *= phi;
        }
    }
}
*/

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
const std::vector<std::vector<Vector>> Mesh::GetInterfaceMidPoints() const { return _cellInterfaceMidPoints; }

double Mesh::GetDistanceToOrigin( unsigned idx_cell ) const {
    double distance = 0.0;
    for( unsigned idx_dim = 0; idx_dim < _dim; idx_dim++ ) {
        distance += _cellMidPoints[idx_cell][idx_dim] * _cellMidPoints[idx_cell][idx_dim];
    }
    return sqrt( distance );
}
