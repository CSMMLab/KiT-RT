#ifndef SETTINGS_H
#define SETTINGS_H

#include <filesystem>
#include <string>

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

    // Mesh
    unsigned _meshDimension;

    // Solver
    unsigned _quadOrder;
    std::string _quadName;
    double _CFL;
    double _tEnd;

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

    // solver Getters
    unsigned GetQuadOrder() const;
    std::string GetQuadName() const;
    double GetCFL() const;
    double GetTEnd() const;

    friend Settings* ReadInputFile( std::string fileName );
};

#endif    // SETTINGS_H
