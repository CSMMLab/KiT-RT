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
                exit( EXIT_FAILURE );
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
                case 2:
                    cellData = vtkDoubleArraySP::New();
                    cellData->SetName( "E(ρU)" );
                    cellData->SetNumberOfComponents( 3 );
                    cellData->SetComponentName( 0, "x" );
                    cellData->SetComponentName( 1, "y" );
                    cellData->SetComponentName( 2, "z" );
                    cellData->SetNumberOfTuples( numCells );
                    for( unsigned j = 0; j < numCells; j++ ) {
                        cellData->SetTuple3( i, results[i][0][j], results[i][1][j], 0.0 );
                    }
                    break;
                case 3:
                    cellData = vtkDoubleArraySP::New();
                    cellData->SetName( "E(ρU)" );
                    cellData->SetNumberOfComponents( 3 );
                    cellData->SetComponentName( 0, "x" );
                    cellData->SetComponentName( 1, "y" );
                    cellData->SetComponentName( 2, "z" );
                    cellData->SetNumberOfTuples( numCells );
                    for( unsigned j = 0; j < numCells; j++ ) {
                        cellData->SetTuple3( i, results[i][0][j], results[i][1][j], results[i][2][j] );
                    }
                    break;
                default:
                    auto log = spdlog::get( "event" );
                    log->error( "[ERROR][IO::ExportVTK] Invalid dimension" );
                    exit( EXIT_FAILURE );
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
    }
    MPI_Barrier( MPI_COMM_WORLD );
}

Mesh* LoadSU2MeshFromFile( const Config* settings ) {
    auto log = spdlog::get( "event" );

    unsigned dim;
    std::vector<Vector> nodes;
    std::vector<std::vector<unsigned>> cells;
    std::vector<std::pair<BOUNDARY_TYPE, std::vector<unsigned>>> boundaries;

    if( !std::filesystem::exists( settings->GetMeshFile() ) ) exit( EXIT_FAILURE );
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
                            exit( EXIT_FAILURE );
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
                            log->error( "[LoadSU2MeshFromFile] Unsupported mesh type!" );
                            exit( EXIT_FAILURE );
                        }
                    }
                }
                break;
            }
        }
        bool mixedElementMesh = !std::equal( numNodesPerCell.begin() + 1, numNodesPerCell.end(), numNodesPerCell.begin() );
        if( mixedElementMesh ) {
            log->error( "[LoadSU2MeshFromFile] Mixed element meshes are currently not supported!" );
            exit( EXIT_FAILURE );
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
        log->error( "[LoadSU2MeshFromFile] File not found" );
        exit( EXIT_FAILURE );
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
                std::cerr << "[ERROR] Unable to open inputfile '" << inputFile << "'!" << std::endl;
                exit( EXIT_FAILURE );
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
