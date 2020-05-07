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

    const unsigned _dim;
    const unsigned _numCells;
    const unsigned _numNodes;
    const unsigned _numNodesPerCell;
    const unsigned _numBoundaries;
    const unsigned _ghostCellID;

    std::vector<Vector> _nodes;
    std::vector<std::vector<unsigned>> _cells;
    std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> _boundaries;

    std::vector<double> _cellAreas;
    std::vector<Vector> _cellMidPoints;
    std::vector<std::vector<unsigned>> _cellNeighbors;
    std::vector<std::vector<Vector>> _cellNormals;
    std::vector<bool> _isBoundaryCell;
    blaze::CompressedMatrix<bool> _nodeNeighbors;
    std::vector<unsigned> _colors;

    void ComputeCellAreas();
    void ComputeCellMidpoints();
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
    const std::vector<Vector>& GetCellMidPoints() const;
    const std::vector<std::vector<unsigned>>& GetCells() const;
    const std::vector<double>& GetCellAreas() const;
    const std::vector<unsigned>& GetPartitionIDs() const;
    const std::vector<std::vector<unsigned>>& GetNeighbours() const;
    const std::vector<std::vector<Vector>>& GetNormals() const;

    /**
     * @brief GetBoundaryCellArray returns array with (dirichlet?) boundary cells
     */
    const std::vector<bool>& GetBoundaryCellArray() const;

    /**
     * @brief ComputeSlopes calculates the slope in every cell into x and y direction
     * @param nq is number of quadrature points
     * @param psiDerX is slope in x direction
     * @param psiDerY is slope in y direction
     * @param psi is solution for which slope is computed
     */
    void ComputeSlopes( unsigned nq, VectorVector& psiDerX, VectorVector& psiDerY, const VectorVector& psi ) const;
};

#endif    // MESH_H
