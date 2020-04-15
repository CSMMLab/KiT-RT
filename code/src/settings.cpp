#include "settings.h"

Settings::Settings() {}

std::string Settings::GetInputFile() const { return _inputFile.string(); }
std::string Settings::GetInputDir() const { return _inputDir.string(); }
std::string Settings::GetOutputFile() const { return _outputFile.string(); }
std::string Settings::GetOutputDir() const { return _outputDir.string(); }
std::string Settings::GetLogDir() const { return _logDir.string(); }
std::string Settings::GetMeshFile() const { return GetInputDir() + _meshFile.string(); }
BOUNDARY_TYPE Settings::GetBoundaryType( std::string name ) const {
    for( unsigned i = 0; i < _boundaries.size(); ++i ) {
        if( name == _boundaries[i].first ) return _boundaries[i].second;
    }
    return BOUNDARY_TYPE::INVALID;
}

unsigned Settings::GetQuadOrder() const { return _quadOrder; }
QUAD_NAME Settings::GetQuadName() const { return _quadName; }
double Settings::GetCFL() const { return _CFL; }
double Settings::GetTEnd() const { return _tEnd; }
