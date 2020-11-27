/*!
 * @file io.h
 * @brief Creates required input from given config file
 */
#ifndef IO_H
#define IO_H

#include "common/typedef.h"
#include <string>
#include <vector>

// Forward Declarations
class Config;
class Mesh;

/*!
 * @brief 
 * @param[in] fileName -name of the file that results shall be exported to
 * @param[in] results - results of computation
 * @param[in] fieldNames - specify different output values
 * @param[in] mesh - mesh results were computed on
 */
void ExportVTK( const std::string fileName,
                const std::vector<std::vector<std::vector<double>>>& results,
                const std::vector<std::string> fieldNames,
                const Mesh* mesh );

/*!
 * @brief Loads mesh from an SU2 file, returns Mesh element
 */
Mesh* LoadSU2MeshFromFile( const Config* settings );


/*!
 * @brief Parses arguments given when calling program from command line
 * @param[in] argc - number arguments
 * @param[in] argv - string with arguments given
 */
std::string ParseArguments( int argc, char* argv[] );

/*!
 * @brief prints configurations set in inputFile to console
 */
void PrintLogHeader( std::string inputFile );

/*!
 * @brief creates a mesh from image given in imageName and writes it to SU2Filename
 * @param[in] imageName - name of image file
 * @param[in] SU2Filename - name of SU2 file where output is saved
 */
Matrix createSU2MeshFromImage( std::string imageName, std::string SU2Filename );

#endif    // IO_H
