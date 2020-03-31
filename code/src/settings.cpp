#include "settings.h"

Settings::Settings() {}

std::string Settings::GetInputFile() const { return _inputFile.string(); }
std::string Settings::GetInputDir() const { return _inputDir.string(); }
std::string Settings::GetOutputFile() const { return _outputFile.string(); }
std::string Settings::GetOutputDir() const { return _outputDir.string(); }
std::string Settings::GetLogDir() const { return _logDir.string(); }

unsigned Settings::GetQuadOrder() const { return _quadOrder; }
std::string Settings::GetQuadName() const { return _quadName; }
