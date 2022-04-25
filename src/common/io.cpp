/*!
 * \file io.cpp
 * \brief Set of utility io functions for rtns
 * \author  J. Kusch, S. Schotthoefer, P. Stammer,  J. Wolters, T. Xiao
 */

#include "common/io.hpp"
#include "common/config.hpp"
#include "common/mesh.hpp"
#include "common/typedef.hpp"
#include "toolboxes/errormessages.hpp"
#include "toolboxes/textprocessingtoolbox.hpp"

#include <iostream>

#include <mpi.h>
#include <omp.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkQuad.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

using vtkPointsSP                 = vtkSmartPointer<vtkPoints>;
using vtkUnstructuredGridSP       = vtkSmartPointer<vtkUnstructuredGrid>;
using vtkTriangleSP               = vtkSmartPointer<vtkTriangle>;
using vtkCellArraySP              = vtkSmartPointer<vtkCellArray>;
using vtkDoubleArraySP            = vtkSmartPointer<vtkDoubleArray>;
using vtkUnstructuredGridWriterSP = vtkSmartPointer<vtkUnstructuredGridWriter>;
using vtkCellDataToPointDataSP    = vtkSmartPointer<vtkCellDataToPointData>;

void ExportVTK( const std::string fileName,
                const std::vector<std::vector<std::vector<double>>>& outputFields,
                const std::vector<std::vector<std::string>>& outputFieldNames,
                const Mesh* mesh ) {
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) {
        unsigned dim             = mesh->GetDim();
        unsigned numCells        = mesh->GetNumCells();
        unsigned numNodes        = mesh->GetNumNodes();
        auto nodes               = mesh->GetNodes();
        auto cells               = mesh->GetCells();
        unsigned numNodesPerCell = mesh->GetNumNodesPerCell();

        auto writer                 = vtkUnstructuredGridWriterSP::New();
        std::string fileNameWithExt = fileName;
        if( !TextProcessingToolbox::StringEndsWith( fileNameWithExt, ".vtk" ) ) {
            fileNameWithExt.append( ".vtk" );
        }
        writer->SetFileName( fileNameWithExt.c_str() );
        auto grid = vtkUnstructuredGridSP::New();
        auto pts  = vtkPointsSP::New();
        pts->SetDataTypeToDouble();
        pts->SetNumberOfPoints( static_cast<int>( numNodes ) );
        unsigned nodeID = 0;
        for( const auto& node : nodes ) {
            if( dim == 2 ) {
                pts->SetPoint( nodeID++, node[0], node[1], 0.0 );
            }
            else if( dim == 3 ) {
                pts->SetPoint( nodeID++, node[0], node[1], node[2] );
            }
            else {
                ErrorMessages::Error( "Unsupported dimension (d=" + std::to_string( dim ) + ")!", CURRENT_FUNCTION );
            }
        }
        vtkCellArraySP cellArray = vtkCellArraySP::New();
        for( unsigned i = 0; i < numCells; ++i ) {
            if( numNodesPerCell == 3 ) {
                auto tri = vtkTriangleSP::New();
                for( unsigned j = 0; j < numNodesPerCell; ++j ) {
                    tri->GetPointIds()->SetId( j, cells[i][j] );
                }
                cellArray->InsertNextCell( tri );
            }
            if( numNodesPerCell == 4 ) {
                auto quad = vtkQuad::New();
                for( unsigned j = 0; j < numNodesPerCell; ++j ) {
                    quad->GetPointIds()->SetId( j, cells[i][j] );
                }
                cellArray->InsertNextCell( quad );
            }
        }
        if( numNodesPerCell == 3 ) {
            grid->SetCells( VTK_TRIANGLE, cellArray );
        }
        else if( numNodesPerCell == 4 ) {
            grid->SetCells( VTK_QUAD, cellArray );
        }

        // Write the output
        for( unsigned idx_group = 0; idx_group < outputFields.size(); idx_group++ ) {

            for( unsigned idx_field = 0; idx_field < outputFields[idx_group].size(); idx_field++ ) {    // Loop over all output fields

                auto cellData = vtkDoubleArraySP::New();
                cellData->SetName( outputFieldNames[idx_group][idx_field].c_str() );

                for( unsigned idx_cell = 0; idx_cell < numCells; idx_cell++ ) {
                    cellData->InsertNextValue( outputFields[idx_group][idx_field][idx_cell] );
                }

                grid->GetCellData()->AddArray( cellData );
            }
        }

        grid->SetPoints( pts );
        grid->Squeeze();

        auto converter = vtkCellDataToPointDataSP::New();
        converter->AddInputDataObject( grid );
        converter->PassCellDataOn();
        converter->Update();

        auto conv_grid = converter->GetOutput();

        writer->SetInputData( conv_grid );

        writer->Write();

        // auto log = spdlog::get( "event" );
        // log->info( "Result successfully exported to '{0}'!", fileNameWithExt );
    }
    MPI_Barrier( MPI_COMM_WORLD );
}

Mesh* LoadSU2MeshFromFile( const Config* settings ) {
    auto log = spdlog::get( "event" );
    log->info( "| Importing mesh. This may take a while for large meshes." );
    unsigned dim;
    std::vector<Vector> nodes;
    std::vector<std::vector<unsigned>> cells;
    std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> boundaries;

    if( !std::filesystem::exists( settings->GetMeshFile() ) )
        ErrorMessages::Error( "Cannot find mesh file '" + settings->GetMeshFile() + "!", CURRENT_FUNCTION );
    std::ifstream ifs( settings->GetMeshFile(), std::ios::in );
    std::string line;

    if( ifs.is_open() ) {
        while( getline( ifs, line ) ) {
            if( line.find( "NDIME", 0 ) != std::string::npos ) {
                dim = static_cast<unsigned>( TextProcessingToolbox::GetTrailingNumber( line ) );
                // if( settings->GetDim() != dim ) {
                //     log->info( "Warning: Mesh dimension does not coinside with problem dimension! Proceed with caution!" );
                // }
                break;
            }
        }
        ifs.clear();
        ifs.seekg( 0, std::ios::beg );
        while( getline( ifs, line ) ) {
            if( line.find( "NPOIN", 0 ) != std::string::npos ) {
                unsigned numPoints = static_cast<unsigned>( TextProcessingToolbox::GetTrailingNumber( line ) );
                nodes.resize( numPoints, Vector( dim, 0.0 ) );
                for( unsigned i = 0; i < numPoints; ++i ) {
                    getline( ifs, line );
                    std::stringstream ss;
                    ss << line;
                    unsigned id = 0;
                    Vector coords( dim, 0.0 );
                    while( !ss.eof() ) {
                        for( unsigned d = 0; d < dim; ++d ) {
                            ss >> coords[d];
                        }
                        ss >> id;
                    }

                    nodes[id] = Vector( coords );
                }
                break;
            }
        }
        ifs.clear();
        ifs.seekg( 0, std::ios::beg );
        while( getline( ifs, line ) ) {
            if( line.find( "NMARK", 0 ) != std::string::npos ) {
                unsigned numBCs = static_cast<unsigned>( TextProcessingToolbox::GetTrailingNumber( line ) );
                boundaries.resize( numBCs );
                for( unsigned i = 0; i < numBCs; ++i ) {
                    std::string markerTag;
                    BOUNDARY_TYPE btype;
                    std::vector<unsigned> bnodes;
                    for( unsigned k = 0; k < 2; ++k ) {
                        getline( ifs, line );
                        if( line.find( "MARKER_TAG", 0 ) != std::string::npos ) {
                            markerTag    = line.substr( line.find( "=" ) + 1 );
                            auto end_pos = std::remove_if(
                                markerTag.begin(), markerTag.end(), []( char c ) { return std::isspace( static_cast<unsigned char>( c ) ); } );
                            markerTag.erase( end_pos, markerTag.end() );
                            btype = settings->GetBoundaryType( markerTag );
                            if( btype == BOUNDARY_TYPE::INVALID ) {
                                std::string errorMsg = std::string( "Invalid Boundary at marker \"" + markerTag + "\"." );
                                ErrorMessages::Error( errorMsg, CURRENT_FUNCTION );
                            }
                        }
                        else if( line.find( "MARKER_ELEMS", 0 ) != std::string::npos ) {
                            unsigned numMarkerElements = static_cast<unsigned>( TextProcessingToolbox::GetTrailingNumber( line ) );
                            for( unsigned j = 0; j < numMarkerElements; ++j ) {
                                getline( ifs, line );
                                std::stringstream ss;
                                ss << line;
                                unsigned type = 0, id = 0;
                                while( !ss.eof() ) {
                                    ss >> type;
                                    for( unsigned d = 0; d < dim; ++d ) {
                                        ss >> id;
                                        bnodes.push_back( id );
                                    }
                                }
                            }
                        }
                        else {
                            ErrorMessages::Error( "Invalid mesh file detected! Make sure boundaries are provided.'", CURRENT_FUNCTION );
                        }
                    }
                    std::sort( bnodes.begin(), bnodes.end() );
                    auto last = std::unique( bnodes.begin(), bnodes.end() );
                    bnodes.erase( last, bnodes.end() );
                    boundaries[i] = std::make_pair( btype, bnodes );
                }
                break;
            }
        }
        ifs.clear();
        ifs.seekg( 0, std::ios::beg );
        std::vector<unsigned> numNodesPerCell;
        while( getline( ifs, line ) ) {
            if( line.find( "NELEM", 0 ) != std::string::npos ) {
                unsigned numCells = static_cast<unsigned>( TextProcessingToolbox::GetTrailingNumber( line ) );
                numNodesPerCell.resize( numCells, 0u );
                for( unsigned i = 0; i < numCells; ++i ) {
                    getline( ifs, line );
                    std::stringstream ss;
                    ss << line;
                    unsigned type = 0;
                    ss >> type;
                    switch( type ) {
                        case 5: {
                            numNodesPerCell[i] = 3;
                            break;
                        }
                        case 9: {
                            numNodesPerCell[i] = 4;
                            break;
                        }
                        default: {
                            ErrorMessages::Error( "Unsupported mesh type!'", CURRENT_FUNCTION );
                        }
                    }
                }
                break;
            }
        }
        bool mixedElementMesh = !std::equal( numNodesPerCell.begin() + 1, numNodesPerCell.end(), numNodesPerCell.begin() );
        if( mixedElementMesh ) {
            ErrorMessages::Error( "Mixed element meshes are currently not supported!'", CURRENT_FUNCTION );
        }
        ifs.clear();
        ifs.seekg( 0, std::ios::beg );
        while( getline( ifs, line ) ) {
            if( line.find( "NELEM", 0 ) != std::string::npos ) {
                unsigned numCells = static_cast<unsigned>( TextProcessingToolbox::GetTrailingNumber( line ) );
                cells.resize( numCells, std::vector<unsigned>( numNodesPerCell[0], 0u ) );
                for( unsigned i = 0; i < numCells; ++i ) {
                    getline( ifs, line );
                    std::stringstream ss;
                    ss << line;
                    unsigned type = 0, id = 0;
                    std::vector<unsigned> buffer( numNodesPerCell[0], 0u );
                    while( !ss.eof() ) {
                        ss >> type;
                        for( unsigned d = 0; d < numNodesPerCell[i]; ++d ) {
                            ss >> buffer[d];
                        }
                        ss >> id;
                    }
                    std::copy( buffer.begin(), buffer.end(), cells[id].begin() );
                }
                break;
            }
        }
    }
    else {
        ErrorMessages::Error( "Cannot open mesh file '" + settings->GetMeshFile() + "!", CURRENT_FUNCTION );
    }
    ifs.close();
    // log->info( "| Mesh imported." );
    return new Mesh( nodes, cells, boundaries );
}

std::string ParseArguments( int argc, char* argv[] ) {
    std::string inputFile;
    std::string usage_help = "\n"
                             "Usage: " +
                             std::string( argv[0] ) + " inputfile\n";

    if( argc < 2 ) {
        std::cout << usage_help;
        exit( EXIT_FAILURE );
    }
    for( int i = 1; i < argc; i++ ) {
        std::string arg = argv[i];
        if( arg == "-h" ) {
            std::cout << usage_help;
            exit( EXIT_SUCCESS );
        }
        else {
            inputFile = std::string( argv[i] );
            std::ifstream f( inputFile );
            if( !f.is_open() ) {
                ErrorMessages::Error( "Unable to open inputfile '" + inputFile + "' !", CURRENT_FUNCTION );
            }
        }
    }
    return inputFile;
}

void PrintLogHeader( std::string inputFile ) {
    int nprocs, rank;
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) {
        auto log = spdlog::get( "event" );

        // New design
        log->info( "------------------------------------------------------------------------" );
        log->info( "|    _  _____ _____     ____ _____                                     |" );
        log->info( "|   | |/ /_ _|_   _|   |  _ \\_   _|                                    |" );
        log->info( "|   | ' / | |  | |_____| |_) || |            Version                   |" );
        log->info( "|   | . \\ | |  | |_____|  _ < | |             1.0.0                    |" );
        log->info( "|   |_|\\_\\___| |_|     |_| \\_\\|_|                                      |" );
        log->info( "|                                                                      |" );
        log->info( "------------------------------------------------------------------------" );
        log->info( "|  Copyright: MIT Licence                                              |" );
        // log->info( "|  Authors: S. Schotthöfer, J. Kusch, J. Wolters, P. Stammer ad T. Xiao |" );
        log->info( "|  Date: 21th April 2022                                               |" );
        log->info( "------------------------------------------------------------------------" );
        log->info( "|" );
        log->info( "| Git commit :\t{0}", GIT_HASH );
        log->info( "| Config file:\t{0}", inputFile );
        log->info( "| MPI Threads:\t{0}", nprocs );
        log->info( "| OMP Threads:\t{0}", omp_get_max_threads() );
        log->info( "|" );
        log->info( "-------------------------- Config File Info ----------------------------" );
        log->info( "|" );
        // print file content while omitting comments
        std::ifstream ifs( inputFile );
        if( ifs.is_open() ) {
            std::string line;
            while( !ifs.eof() ) {
                std::getline( ifs, line );
                if( line[0] != '%' ) log->info( "| {0}", line );
            }
        }
        // log->info( "------------------------------------------------------------------------" );
    }
    MPI_Barrier( MPI_COMM_WORLD );
}

Matrix createSU2MeshFromImage( std::string imageName, std::string SU2Filename ) {
    auto log = spdlog::get( "event" );

    // std::cout << "imageName" << imageName << "\n";
    // std::cout << "pyhtonPath" << KITRT_PYTHON_PATH;

    if( !std::filesystem::exists( imageName ) ) {
        ErrorMessages::Error( "Can not open image '" + imageName + "'!", CURRENT_FUNCTION );
    }
    std::filesystem::path outDir( std::filesystem::path( SU2Filename ).parent_path() );
    if( !std::filesystem::exists( outDir ) ) {
        ErrorMessages::Error( "Output directory '" + outDir.string() + "' does not exists!", CURRENT_FUNCTION );
    }

    std::string pyPath = KITRT_PYTHON_PATH;

    if( !Py_IsInitialized() ) {
        Py_InitializeEx( 0 );
        if( !Py_IsInitialized() ) {
            ErrorMessages::Error( "Python init failed!", CURRENT_FUNCTION );
        }
        PyRun_SimpleString( ( "import sys\nsys.path.append('" + pyPath + "')" ).c_str() );
    }

    PyObject *pArgs, *pReturn, *pModule, *pFunc;
    PyArrayObject* np_ret;

    auto image = PyUnicode_FromString( imageName.c_str() );
    auto su2   = PyUnicode_FromString( SU2Filename.c_str() );

    std::string moduleName = "mesh_from_image";
    pModule                = PyImport_ImportModule( moduleName.c_str() );
    if( !pModule ) {
        PyErr_Print();
        ErrorMessages::Error( "'" + moduleName + "' can not be imported!", CURRENT_FUNCTION );
    }

    pFunc = PyObject_GetAttrString( pModule, "generate" );
    if( !pFunc || !PyCallable_Check( pFunc ) ) {
        PyErr_Print();
        Py_DecRef( pModule );
        Py_DecRef( pFunc );
        ErrorMessages::Error( "'generate' is null or not callable!", CURRENT_FUNCTION );
    }

    pArgs = PyTuple_New( 2 );
    PyTuple_SetItem( pArgs, 0, reinterpret_cast<PyObject*>( image ) );
    PyTuple_SetItem( pArgs, 1, reinterpret_cast<PyObject*>( su2 ) );
    pReturn = PyObject_CallObject( pFunc, pArgs );
    np_ret  = reinterpret_cast<PyArrayObject*>( pReturn );

    size_t m{ static_cast<size_t>( PyArray_SHAPE( np_ret )[0] ) };
    size_t n{ static_cast<size_t>( PyArray_SHAPE( np_ret )[1] ) };
    double* c_out = reinterpret_cast<double*>( PyArray_DATA( np_ret ) );

    Matrix gsImage( m, n, c_out );

    // Finalizing
    Py_DecRef( pFunc );
    Py_DecRef( pModule );
    Py_DECREF( np_ret );

    return gsImage.transpose();
}
