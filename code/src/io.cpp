#include "io.h"
#include "toolboxes/errormessages.h"

void ExportVTK( const std::string fileName,
                const std::vector<std::vector<std::vector<double>>>& results,
                const std::vector<std::string> fieldNames,
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
        if( fileNameWithExt.substr( fileNameWithExt.find_last_of( "." ) + 1 ) != ".vtk" ) {
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

        for( unsigned i = 0; i < results.size(); i++ ) {
            auto cellData = vtkDoubleArraySP::New();
            cellData->SetName( fieldNames[i].c_str() );
            switch( results[i].size() ) {
                case 1:
                    for( unsigned j = 0; j < numCells; j++ ) {
                        cellData->InsertNextValue( results[i][0][j] );
                    }
                    break;
                default:
                    auto log = spdlog::get( "event" );
                    ErrorMessages::Error( "Please implement output for results of size " + std::to_string( results[i].size() ) + "!",
                                          CURRENT_FUNCTION );
            }
            grid->GetCellData()->AddArray( cellData );
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

        auto log = spdlog::get( "event" );
        log->info( "Result successfully exported to '{0}'!", fileNameWithExt );
    }
    MPI_Barrier( MPI_COMM_WORLD );
}

Mesh* LoadSU2MeshFromFile( const Config* settings ) {
    auto log = spdlog::get( "event" );

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
                dim = static_cast<unsigned>( std::stoi( line.substr( line.find_first_of( "0123456789" ), line.length() - 1 ) ) );
                break;
            }
        }
        ifs.clear();
        ifs.seekg( 0, std::ios::beg );
        while( getline( ifs, line ) ) {
            if( line.find( "NPOIN", 0 ) != std::string::npos ) {
                unsigned numPoints = static_cast<unsigned>( std::stoi( line.substr( line.find_first_of( "0123456789" ), line.length() - 1 ) ) );
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
                unsigned numBCs = static_cast<unsigned>( std::stoi( line.substr( line.find_first_of( "0123456789" ), line.length() - 1 ) ) );
                boundaries.resize( numBCs );
                for( unsigned i = 0; i < numBCs; ++i ) {
                    std::string markerTag;
                    BOUNDARY_TYPE btype;
                    std::vector<unsigned> bnodes;
                    for( unsigned k = 0; k < 2; ++k ) {
                        getline( ifs, line );
                        if( line.find( "MARKER_TAG", 0 ) != std::string::npos ) {
                            markerTag    = line.substr( line.find( "=" ) + 1 );
                            auto end_pos = std::remove_if( markerTag.begin(), markerTag.end(), isspace );
                            markerTag.erase( end_pos, markerTag.end() );
                            btype = settings->GetBoundaryType( markerTag );
                            if( btype == BOUNDARY_TYPE::INVALID ) {
                                std::string errorMsg = std::string( "Invalid Boundary at marker \"" + markerTag + "\"." );
                                ErrorMessages::Error( errorMsg, CURRENT_FUNCTION );
                            }
                        }
                        else if( line.find( "MARKER_ELEMS", 0 ) != std::string::npos ) {
                            unsigned numMarkerElements =
                                static_cast<unsigned>( std::stoi( line.substr( line.find_first_of( "0123456789" ), line.length() - 1 ) ) );
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
                unsigned numCells = static_cast<unsigned>( std::stoi( line.substr( line.find_first_of( "0123456789" ), line.length() - 1 ) ) );
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
                unsigned numCells = static_cast<unsigned>( std::stoi( line.substr( line.find_first_of( "0123456789" ), line.length() - 1 ) ) );
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
                ErrorMessages::OptionNotSetError( "Unable to open inputfile '" + inputFile + "' !", CURRENT_FUNCTION );
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

        log->info( "RTSN" );
        log->info( "================================================================" );
        log->info( "Git commit :\t{0}", GIT_HASH );
        log->info( "Config file:\t{0}", inputFile );
        log->info( "MPI Threads:\t{0}", nprocs );
        log->info( "OMP Threads:\t{0}", omp_get_max_threads() );
        log->info( "================================================================" );
        // print file content while omitting comments
        std::ifstream ifs( inputFile );
        if( ifs.is_open() ) {
            std::string line;
            while( !ifs.eof() ) {
                std::getline( ifs, line );
                if( line[0] != '%' ) log->info( " {0}", line );
            }
        }
        log->info( "================================================================" );
        log->info( "" );
    }
    MPI_Barrier( MPI_COMM_WORLD );
}

Matrix createSU2MeshFromImage( std::string imageName, std::string SU2Filename ) {
    auto log = spdlog::get( "event" );

    if( !std::filesystem::exists( imageName ) ) {
        ErrorMessages::Error( "Can not open image '" + imageName + "'!", CURRENT_FUNCTION );
        exit( EXIT_FAILURE );
    }
    std::filesystem::path outDir( std::filesystem::path( SU2Filename ).parent_path() );
    if( !std::filesystem::exists( outDir ) ) {
        ErrorMessages::Error( "Output directory '" + outDir.string() + "' does not exists!", CURRENT_FUNCTION );
        exit( EXIT_FAILURE );
    }

    std::cout << getenv( "PYTHONPATH" ) << std::endl;
    setenv( "PYTHONPATH", RTSN_PYTHON_PATH, 1 );
    std::cout << getenv( "PYTHONPATH" ) << std::endl << std::flush;
    Py_Initialize();
    PyObject *pArgs, *pReturn, *pModule, *pFunc;
    PyArrayObject* np_ret;

    auto image = PyUnicode_FromString( imageName.c_str() );
    auto su2   = PyUnicode_FromString( SU2Filename.c_str() );

    PyObject* pName = PyUnicode_FromString( "mesh_from_image" );
    pModule         = PyImport_Import( pName );
    Py_DECREF( pName );
    if( !pModule ) {
        ErrorMessages::Error( "'mesh_from_image.py' can not be imported!", CURRENT_FUNCTION );
        exit( EXIT_FAILURE );
    }

    pFunc = PyObject_GetAttrString( pModule, "generate" );
    if( !pFunc || !PyCallable_Check( pFunc ) ) {
        Py_DECREF( pModule );
        Py_XDECREF( pFunc );
        ErrorMessages::Error( "'generate' is null or not callable!", CURRENT_FUNCTION );
        exit( EXIT_FAILURE );
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
    Py_DECREF( pFunc );
    Py_DECREF( pModule );
    Py_DECREF( np_ret );
    Py_Finalize();

    return gsImage;
}
