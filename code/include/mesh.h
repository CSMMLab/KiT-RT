#ifndef MESH_H
#define MESH_H

#include <vector>

class Mesh
{
  protected:
    unsigned _dim;
    unsigned _numCells;
    unsigned _numNodes;
    unsigned _numNodesPerCell;

    std::vector<std::vector<double>> _nodes;
    std::vector<std::vector<double>> _cells;

  public:
    unsigned GetDim() const;
    unsigned GetNumCells() const;
    unsigned GetNumNodes() const;
    unsigned GetNumNodesPerCell() const;
    const std::vector<std::vector<double>>& GetNodes() const;
    const std::vector<std::vector<double>>& GetCells() const;
};

#endif    // MESH_H
