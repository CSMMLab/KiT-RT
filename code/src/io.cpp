#include "io.h"
#include "toolboxes/CRTSNError.h"

void ExportVTK( const std::string fileName,
                const std::vector<std::vector<std::vector<double>>>& results,
                const std::vector<std::string> fieldNames,
                const CConfig* settings,
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
        writer->SetFileName( ( settings->GetOutputDir() + fileNameWithExt ).c_str() );
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
                default: std::cout << "[ERROR][IO::ExportVTK] Invalid dimension" << std::endl;
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

void InitLogger( std::string logDir, spdlog::level::level_enum terminalLogLvl, spdlog::level::level_enum fileLogLvl ) {
    // create log dir if not existent
    if( !std::filesystem::exists( logDir ) ) {
        std::filesystem::create_directory( logDir );
    }

    // create sinks if level is not off
    std::vector<spdlog::sink_ptr> sinks;
    if( terminalLogLvl != spdlog::level::off ) {
        // create spdlog terminal sink
        auto terminalSink = std::make_shared<spdlog::sinks::stdout_sink_mt>();
        terminalSink->set_level( terminalLogLvl );
        terminalSink->set_pattern( "%v" );
        sinks.push_back( terminalSink );
    }
    if( fileLogLvl != spdlog::level::off ) {
        // define filename on root
        int pe;
        MPI_Comm_rank( MPI_COMM_WORLD, &pe );
        char cfilename[1024];
        if( pe == 0 ) {
            // get date and time
            time_t now = time( nullptr );
            struct tm tstruct;
            char buf[80];
            tstruct = *localtime( &now );
            strftime( buf, sizeof( buf ), "%Y-%m-%d_%X", &tstruct );

            // set filename to date and time
            std::string filename = buf;

            // in case of existing files append '_#'
            int ctr = 0;
            if( std::filesystem::exists( logDir + filename ) ) {
                filename += "_" + std::to_string( ++ctr );
            }
            while( std::filesystem::exists( logDir + filename ) ) {
                filename.pop_back();
                filename += std::to_string( ++ctr );
            }
            strncpy( cfilename, filename.c_str(), sizeof( cfilename ) );
            cfilename[sizeof( cfilename ) - 1] = 0;
        }
        MPI_Bcast( &cfilename, sizeof( cfilename ), MPI_CHAR, 0, MPI_COMM_WORLD );
        MPI_Barrier( MPI_COMM_WORLD );

        // create spdlog file sink
        auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( logDir + cfilename );
        fileSink->set_level( fileLogLvl );
        fileSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );
        sinks.push_back( fileSink );
    }

    // register all sinks
    auto event_logger = std::make_shared<spdlog::logger>( "event", begin( sinks ), end( sinks ) );
    spdlog::register_logger( event_logger );
    spdlog::flush_every( std::chrono::seconds( 5 ) );
}

Mesh* LoadSU2MeshFromFile( const CConfig* settings ) {
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
                                CRTSNError::Error( errorMsg, CURRENT_FUNCTION );
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
                if( line[0] != '#' ) log->info( " {0}", line );
            }
        }
        log->info( "================================================================" );
        log->info( "" );
    }
    MPI_Barrier( MPI_COMM_WORLD );
}

Settings* ReadInputFile( std::string inputFile ) {
    bool validConfig = true;

    Settings* settings = new Settings;

    MPI_Comm_rank( MPI_COMM_WORLD, &settings->_comm_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &settings->_comm_size );

    try {
        auto file = cpptoml::parse_file( inputFile );

        settings->_inputFile = std::filesystem::path( inputFile );

        auto cwd        = std::filesystem::current_path();
        std::string tmp = std::filesystem::path( inputFile ).parent_path().string();
        if( tmp.substr( tmp.size() - 1 ) != "/" ) tmp.append( "/" );
        settings->_inputDir = tmp;

        // section IO
        auto io       = file->get_table( "io" );
        auto meshFile = io->get_as<std::string>( "meshFile" );
        if( meshFile ) {
            settings->_meshFile = std::filesystem::path( *meshFile );
        }
        else {
            spdlog::error( "[inputfile] [io] 'meshFile' not set!" );
            validConfig = false;
        }

        auto outputDir = io->get_as<std::string>( "outputDir" );
        if( outputDir ) {
            std::string tmp = *outputDir;
            if( tmp.substr( tmp.size() - 1 ) != "/" ) tmp.append( "/" );
            settings->_outputDir = std::filesystem::path( tmp );
        }
        else {
            spdlog::error( "[inputfile] [io] 'outputDir' not set!" );
            validConfig = false;
        }

        auto outputFile = io->get_as<std::string>( "outputFile" );
        if( outputFile ) {
            settings->_outputFile = std::filesystem::path( settings->_inputDir.string() + settings->_outputDir.string() + *outputFile );
        }
        else {
            spdlog::error( "[inputfile] [io] 'outputFile' not set!" );
            validConfig = false;
        }

        auto logDir = io->get_as<std::string>( "logDir" );
        if( logDir ) {
            std::string tmp = *logDir;
            if( tmp.substr( tmp.size() - 1 ) != "/" ) tmp.append( "/" );
            settings->_logDir = std::filesystem::path( tmp );
        }
        else {
            spdlog::error( "[inputfile] [io] 'logDir' not set!" );
            validConfig = false;
        }

        // section solver
        auto solver = file->get_table( "solver" );
        auto CFL    = solver->get_as<double>( "CFL" );
        if( CFL ) {
            settings->_CFL = *CFL;
        }
        else {
            spdlog::error( "[inputfile] [solver] 'CFL' not set!" );
            validConfig = false;
        }

        auto tEnd = solver->get_as<double>( "tEnd" );
        if( tEnd ) {
            settings->_tEnd = *tEnd;
        }
        else {
            spdlog::error( "[inputfile] [solver] 'tEnd' not set!" );
            validConfig = false;
        }

        auto quadType = solver->get_as<std::string>( "quadType" );
        if( quadType ) {
            std::string quadTypeString = *quadType;
            try {
                settings->_quadName = Quadrature_Map.at( quadTypeString );
            } catch( const std::exception& e ) {
                spdlog::error( "Error: '{0}' is not a feasible quadrature type. Please check the config template!", quadTypeString );
                exit( EXIT_FAILURE );    // Quit RTSN
            }
        }
        else {
            spdlog::error( "[inputfile] [solver] 'quadType' not set!" );
            validConfig = false;
        }

        auto quadOrder = solver->get_as<unsigned>( "quadOrder" );
        if( quadOrder ) {
            settings->_quadOrder = *quadOrder;
        }
        else {
            spdlog::error( "[inputfile] [solver] 'quadOrder' not set!" );
            validConfig = false;
        }

        auto BCStrings = solver->get_array_of<cpptoml::array>( "boundaryConditions" );
        if( BCStrings ) {
            for( unsigned i = 0; i < BCStrings->size(); ++i ) {
                auto BCString = ( *BCStrings )[i]->get_array_of<std::string>();
                BOUNDARY_TYPE type;
                if( ( *BCString )[1].compare( "dirichlet" ) == 0 ) {
                    type = BOUNDARY_TYPE::DIRICHLET;
                }
                else {
                    spdlog::error( "[inputfile] [solver] Encountered invalid boundary type '" + ( *BCString )[0] + "'!" );
                    exit( EXIT_FAILURE );
                }
                settings->_boundaries.push_back( std::make_pair( ( *BCString )[0], type ) );
            }
        }
        else {
            spdlog::error( "[inputfile] [solver] 'boundaryConditions' Not set!" );
            exit( EXIT_FAILURE );
        }

    } catch( const cpptoml::parse_exception& e ) {
        spdlog::error( "Failed to parse {0}: {1}", inputFile.c_str(), e.what() );
        exit( EXIT_FAILURE );
    }

    if( !validConfig ) {
        exit( EXIT_FAILURE );
    }

    return settings;
}
