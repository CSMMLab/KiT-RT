#ifndef IO_H
#define IO_H

#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>

#include "cpptoml.h"
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
#include "settings/config.h"

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
                const Mesh* mesh );

Mesh* LoadSU2MeshFromFile( const Config* settings );

std::string ParseArguments( int argc, char* argv[] );

void PrintLogHeader( std::string inputFile );

#endif    // IO_H
