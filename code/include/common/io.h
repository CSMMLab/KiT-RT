#ifndef IO_H
#define IO_H

#include "common/typedef.h"
#include <string>
#include <vector>

// Forward Declarations
class Config;
class Mesh;
// class Matrix;

void ExportVTK( const std::string fileName,
                const std::vector<std::vector<std::vector<double>>>& results,
                const std::vector<std::string> fieldNames,
                const Mesh* mesh );

Mesh* LoadSU2MeshFromFile( const Config* settings );

std::string ParseArguments( int argc, char* argv[] );

void PrintLogHeader( std::string inputFile );

Matrix createSU2MeshFromImage( std::string imageName, std::string SU2Filename );

#endif    // IO_H
