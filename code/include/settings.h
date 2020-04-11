#ifndef SETTINGS_H
#define SETTINGS_H

#include <filesystem>
#include <string>

#include "mesh.h"

class Settings
{
  private:
    // IO
    std::filesystem::path _inputFile;
    std::filesystem::path _inputDir;
    std::filesystem::path _outputFile;
    std::filesystem::path _outputDir;
    std::filesystem::path _logDir;
    std::filesystem::path _meshFile;

    std::vector<std::pair<std::string, BOUNDARY_TYPE>> _boundaries;

    // Mesh
    unsigned _meshDimension;

    // Solver
    unsigned _quadOrder;
    double _CFL;

    // MPI
    int _comm_rank;
    int _comm_size;

  public:
    Settings();

    std::string GetInputFile() const;
    std::string GetInputDir() const;
    std::string GetOutputFile() const;
    std::string GetOutputDir() const;
    std::string GetLogDir() const;
    std::string GetMeshFile() const;

    BOUNDARY_TYPE GetBoundaryType( std::string name ) const;

    friend Settings* ReadInputFile( std::string fileName );
};

#endif    // SETTINGS_H
