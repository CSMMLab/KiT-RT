#include "mesh.h"

unsigned Mesh::GetDim() const { return _dim; }

unsigned Mesh::GetNumCells() const { return _numCells; }

unsigned Mesh::GetNumNodes() const { return _numNodes; }

unsigned Mesh::GetNumNodesPerCell() const { return _numNodesPerCell; }

const std::vector<std::vector<double>>& Mesh::GetNodes() const { return _nodes; }

const std::vector<std::vector<double>>& Mesh::GetCells() const { return _cells; }
