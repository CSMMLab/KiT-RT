#ifndef MESH_H
#define MESH_H

#include <algorithm>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "blaze/Blaze.h"
#include "metis.h"
#include "parmetis.h"
#include "spdlog/spdlog.h"

#include "typedef.h"

enum BOUNDARY_TYPE { DIRICHLET, INVALID };

class Mesh
{
  protected:
    std::shared_ptr<spdlog::logger> _log;

    unsigned _dim;
    unsigned _numCells;
    unsigned _numNodes;
    unsigned _numNodesPerCell;
    unsigned _numBoundaries;
    unsigned _ghostCellID;

    std::vector<Vector> _nodes;
    std::vector<std::vector<unsigned>> _cells;
    std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> _boundaries;

    std::vector<double> _cellAreas;
    std::vector<std::vector<unsigned>> _cellNeighbors;
    std::vector<std::vector<Vector>> _cellNormals;
    std::vector<unsigned> _colors;

    void ComputeCellAreas();
    void ComputeConnectivity();
    void ComputeNormals();
    void ComputePartitioning();
    Vector ComputeOutwardFacingNormal( const Vector& nodeA, const Vector& nodeB, const Vector& cellCenter );

  public:
    Mesh() = delete;
    Mesh( std::vector<Vector> nodes,
          std::vector<std::vector<unsigned>> cells,
          std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> boundaries );
    ~Mesh();

    unsigned GetDim() const;
    unsigned GetNumCells() const;
    unsigned GetNumNodes() const;
    unsigned GetNumNodesPerCell() const;
    const std::vector<Vector>& GetNodes() const;
    const std::vector<std::vector<unsigned>>& GetCells() const;
    const std::vector<double>& GetCellAreas() const;
};

#endif    // MESH_H
