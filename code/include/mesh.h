#ifndef MESH_H
#define MESH_H

#include <algorithm>
#include <mpi.h>
#include <omp.h>
#include <vector>

#include "blaze/math/CompressedMatrix.h"
#include "metis.h"
#include "parmetis.h"
#include "spdlog/spdlog.h"

#include "settings/globalconstants.h"
#include "settings/typedef.h"
#include "toolboxes/errormessages.h"

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
    std::vector<BOUNDARY_TYPE> _cellBoundaryTypes;
    blaze::CompressedMatrix<bool> _nodeNeighbors;
    std::vector<unsigned> _colors;

    void ComputeCellAreas();
    void ComputeCellMidpoints();
    void ComputeConnectivity();
    void ComputePartitioning();
    Vector ComputeOutwardFacingNormal( const Vector& nodeA, const Vector& nodeB, const Vector& cellCenter );

  public:
    Mesh() = delete;
    Mesh( std::vector<Vector> nodes,
          std::vector<std::vector<unsigned>> cells,
          std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> boundaries );
    ~Mesh();

    inline unsigned GetDim() const { return _dim; }
    inline unsigned GetNumCells() const { return _numCells; }
    inline unsigned GetNumNodes() const { return _numNodes; }
    inline unsigned GetNumNodesPerCell() const { return _numNodesPerCell; }

    /**
     * @brief Returns all node coordinates
     * @return dimension: numNodes x dim
     */
    const std::vector<Vector>& GetNodes() const;

    /**
     * @brief  Returns the mid point coordinates of each cell
     * @return dimension: numCells x dim
     */
    const std::vector<Vector>& GetCellMidPoints() const;

    /**
     * @brief Returns all node IDs that construct up each cell
     * @return dimension: numCells x numNodes
     */
    const std::vector<std::vector<unsigned>>& GetCells() const;

    /**
     * @brief Returns the cell area of each cell
     * @return dimension: numCells
     */
    const std::vector<double>& GetCellAreas() const;

    /**
     * @brief Return the color/ID of the mesh partition
     * @return dimension: numCells
     */
    const std::vector<unsigned>& GetPartitionIDs() const;

    /**
     * @brief Returns the neighbor cell IDs for every cell
     * @return dimension: numCells x numNodes
     */
    const std::vector<std::vector<unsigned>>& GetNeighbours() const;

    /**
     * @brief Returns the edge length scaled normal vectors of each cell
     * @return dimension: numCells x numNodes x dim
     */
    const std::vector<std::vector<Vector>>& GetNormals() const;

    /**
     * @brief Returns the boundary enum for each cell. BOUNDARY_TYPE::NONE is the default.
     * @return dimension: numCells
     */
    const std::vector<BOUNDARY_TYPE>& GetBoundaryTypes() const;

    double GetDistanceToOrigin( unsigned idx_cell ) const;
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
