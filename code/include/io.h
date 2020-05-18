#ifndef IO_H
#define IO_H

#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>

#include "cpptoml.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "spdlog/spdlog.h"
#include <mpi.h>
#include <omp.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPointDataToCellData.h>
#include <vtkQuad.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>

#include "mesh.h"
#include "settings.h"

#include "settings/CConfig.h"

using vtkPointsSP                 = vtkSmartPointer<vtkPoints>;
using vtkUnstructuredGridSP       = vtkSmartPointer<vtkUnstructuredGrid>;
using vtkTriangleSP               = vtkSmartPointer<vtkTriangle>;
using vtkCellArraySP              = vtkSmartPointer<vtkCellArray>;
using vtkDoubleArraySP            = vtkSmartPointer<vtkDoubleArray>;
using vtkUnstructuredGridWriterSP = vtkSmartPointer<vtkUnstructuredGridWriter>;
using vtkUnstructuredGridReaderSP = vtkSmartPointer<vtkUnstructuredGridReader>;
using vtkCellDataToPointDataSP    = vtkSmartPointer<vtkCellDataToPointData>;
using vtkPointDataToCellDataSP    = vtkSmartPointer<vtkPointDataToCellData>;

void ExportVTK( const std::string fileName,
                const std::vector<std::vector<std::vector<double>>>& results,
                const std::vector<std::string> fieldNames,
                const CConfig* settings,
                const Mesh* mesh );
void InitLogger( std::string logDir, spdlog::level::level_enum terminalLogLvl, spdlog::level::level_enum fileLogLvl );
Mesh* LoadSU2MeshFromFile( const CConfig* settings );
std::string ParseArguments( int argc, char* argv[] );
Settings* ReadInputFile( std::string inputFile );
void PrintLogHeader( std::string inputFile );

#endif    // IO_H
