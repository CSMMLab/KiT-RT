#ifndef MESH_H
#define MESH_H

#include <utility>    // for pair
#include <vector>     // for vector

#include "blaze/math/sparse/CompressedMatrix.h"    // for CompressedMatrix
#include "common/globalconstants.h"                // for BOUNDARY_TYPE
#include "common/typedef.h"                        // for Vector, VectorVector

#include "toolboxes/errormessages.h"
#include "toolboxes/reconstructor.h"

class Mesh
{
  protected:
    const unsigned _dim;             /*!< @brief spatial dimension of the mesh, i.e. 1D,2D,3D */
    const unsigned _numCells;        /*!< @brief number of cells in the mesh */
    const unsigned _numNodes;        /*!< @brief number of nodes in the mesh (for node centered view)*/
    const unsigned _numNodesPerCell; /*!< @brief number of nodes per cell */
    const unsigned _numBoundaries;   /*!< @brief number of boundary cells in the mesh */
    const unsigned _ghostCellID; /*!< @brief Id of the ghost cell. (we use only one ghost cell). equal to _numCells and therefore has the ID of the
                                    last cell + 1 */
    unsigned _numNodesPerBoundary;
    std::vector<std::pair<double, double>> _bounds;    // ???

    std::vector<Vector> _nodes;                /*!< @brief nodes in the mesh. dimension:_numNodes<_dim> */
    std::vector<std::vector<unsigned>> _cells; /*!< @brief cells in the mesh. dimension:_numCells<_numNodesPerCell>  */

    /*! @brief boundary cells in the mesh. Pair defines boundary type of the boundary nodes of the cell. numBoundaries<(1,numBoundaryNodes)>*/
    std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> _boundaries;
    std::vector<double> _cellAreas;                    /*!< @brief cell areas of the mesh. dimension: numCells*/
    std::vector<Vector> _cellMidPoints;                /*!< @brief cell midpoints of the mesh. dimension: numCells<dim>*/
    std::vector<std::vector<unsigned>> _cellNeighbors; /*!< @brief neighbors of each cell. dimension: numCells<numNodesPerCell>*/

    /*! @brief outward facing normals of each side of each cell. dimension: numCells<numNodesPerCell<dim>>, all
                normals are facing away from the cell center, and scaled with the edge length */
    std::vector<std::vector<Vector>> _cellNormals;
    /*! @brief Tags each cell with its boundary type. None means no boundary. dimension: numCells */
    std::vector<BOUNDARY_TYPE> _cellBoundaryTypes;
    blaze::CompressedMatrix<bool> _nodeNeighbors; /*!< @brief neighborshood relationship of nodes for (par-)metis */

    void ComputeCellAreas();     /*!< @brief Computes only the areas of the mesh cells. Write to _cellAreas. */
    void ComputeCellMidpoints(); /*!< @brief Compute only the midpoints of the cells. Write to _cellMidPoints*/
    void ComputeConnectivity();  /*!< @brief Computes _cellNeighbors and _nodeNeighbors, i.e. neighborship relation in mesh*/

    /*! @brief Computes outward facing normal of two neighboring nodes nodeA and nodeB with common cellCellcenter.
     *          Normals are scaled with their respective edge length
     *  @param nodeA: first node
     *  @param nodeB: neighboring node to nodeA
     *  @param cellCenter: Center of the cell that has nodeA and nodeB as nodes.
     *  @return outward facing normal */
    Vector ComputeOutwardFacingNormal( const Vector& nodeA, const Vector& nodeB, const Vector& cellCenter );
    void ComputeBounds(); /*!< @brief Computes the spatial bounds of a 2D domain. */

  public:
    Mesh() = delete;    //  no default constructor

    /*! @brief Constructor of mesh. Needs nodes, cells, and boundary descriptions as specified above.
     *          See LoadSU2MeshFromFile in io.cpp for setup information*/
    Mesh( std::vector<Vector> nodes,
          std::vector<std::vector<unsigned>> cells,
          std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> boundaries );
    ~Mesh();

    inline unsigned GetDim() const { return _dim; }
    inline unsigned GetNumCells() const { return _numCells; }
    inline unsigned GetNumNodes() const { return _numNodes; }
    inline unsigned GetNumNodesPerCell() const { return _numNodesPerCell; }

    /*! @brief Returns all node coordinates
     *  @return dimension: numNodes x dim */
    const std::vector<Vector>& GetNodes() const;

    /*! @brief  Returns the mid point coordinates of each cell
     *  @return dimension: numCells x dim */
    const std::vector<Vector>& GetCellMidPoints() const;

    /*! @brief Returns all node IDs that construct up each cell
     *  @return dimension: numCells x numNodes */
    const std::vector<std::vector<unsigned>>& GetCells() const;

    /*! @brief Returns the cell area of each cell
     *  @return dimension: numCells */
    const std::vector<double>& GetCellAreas() const;

    /*! @brief Returns the neighbor cell IDs for every cell
     *  @return dimension: numCells x numNodes */
    const std::vector<std::vector<unsigned>>& GetNeighbours() const;

    /*! @brief Returns the edge length scaled normal vectors of each cell
     *  @return dimension: numCells x numNodes x dim */
    const std::vector<std::vector<Vector>>& GetNormals() const;

    /*! @brief Returns the boundary enum for each cell. BOUNDARY_TYPE::NONE is the default.
     *  @return dimension: numCells */
    const std::vector<BOUNDARY_TYPE>& GetBoundaryTypes() const;

    /*! @brief Returns the minimal and maximal coordinates of all nodes for each dimension
     *  @return dimension: dim */
    const std::vector<std::pair<double, double>> GetBounds() const;

    /*! @brief Returns distance of a specified cells center to the coordinate systems origin
     *  @return dimension: scalar */
    double GetDistanceToOrigin( unsigned idx_cell ) const;

    /*! @brief ComputeSlopes calculates the slope in every cell into x and y direction using the divergence theorem.
     *  @param nq is number of quadrature points
     *  @param psiDerX is slope in x direction (gets computed. Slope is stored here)
     *  @param psiDerY is slope in y direction (gets computed. Slope is stored here)
     *  @param psi is solution for which slope is computed */
    // Not used
    void ComputeSlopes( unsigned nq, VectorVector& psiDerX, VectorVector& psiDerY, const VectorVector& psi ) const;

    /*! @brief Structured mesh slope reconstruction with flux limiters.
     *  @param nq is number of quadrature points
     *  @param psiDerX is slope in x direction (gets computed. Slope is stored here)
     *  @param psiDerY is slope in y direction (gets computed. Slope is stored here)
     *  @param psi is solution for which slope is computed */
    void ReconstructSlopesS( unsigned nq, VectorVector& psiDerX, VectorVector& psiDerY, const VectorVector& psi ) const;
    /*! @brief Use gauss theorem and limiters. For unstructured mesh *
     *  @param nq is number of quadrature points
     *  @param psiDerX is slope in x direction (gets computed. Slope is stored here)
     *  @param psiDerY is slope in y direction (gets computed. Slope is stored here)
     *  @param psi is solution for which slope is computed */
    void ReconstructSlopesU( unsigned nq, VectorVector& psiDerX, VectorVector& psiDerY, const VectorVector& psi ) const;
};

#endif    // MESH_H
