#ifndef IO_H
#define IO_H

#include "common/typedef.hpp"
#include <string>
#include <vector>

// Forward Declarations
class Config;
class Mesh;

/*! @brief Function to export solver Volume output to VTK file.
 *  @param fileName Filename of output file
 *  @param outputFields numerical output of the solver. Dimensions: (OutputGroupSize, OutputFieldSize, NumberMeshCells)
 *  @param outputFieldNames names of the outputfields. Dimensions: (OutputGroupSize, OutputFieldSize)
 *  @param mesh Mesh with cells (the mesh used for the computation)
 */
void ExportVTK( const std::string fileName,
                const std::vector<std::vector<std::vector<double>>>& outputFields,
                const std::vector<std::vector<std::string>>& outputFieldNames,
                const Mesh* mesh );

Mesh* LoadSU2MeshFromFile( const Config* settings );

std::string ParseArguments( int argc, char* argv[] );

void PrintLogHeader( std::string inputFile );

//Matrix createSU2MeshFromImage( std::string imageName, std::string SU2Filename ); Deprecated

#endif    // IO_H
