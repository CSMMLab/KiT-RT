#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow() : QMainWindow( nullptr ), ui( new Ui::MainWindow ) { initUI(); }

MainWindow::~MainWindow() { delete ui; }

void MainWindow::initUI() {
    ui->setupUi( this );

    this->statusBar()->setStyleSheet( "color:white; background-color: rgb(40, 40, 40);" );
    this->setStyleSheet( "color:white; background-color: rgb(50, 50, 50);" );

    QHBoxLayout* outputHB = new QHBoxLayout( ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( new QLabel( "Output directory:", ui->settingsFrame ) );
    _outputDir = new QLineEdit( ui->settingsFrame );
    _outputDir->setPlaceholderText( "/path/to/dir" );
    outputHB->addWidget( _outputDir );
    QPushButton* outputDirButton = new QPushButton( "Browse", ui->settingsFrame );
    outputHB->addWidget( outputDirButton );
    ui->settingsFrame->layout()->addItem( outputHB );
    connect( outputDirButton, SIGNAL( clicked() ), this, SLOT( selectOutputDir() ) );

    ui->settingsFrame->layout()->addWidget( new QLabel( "Output filename:", ui->settingsFrame ) );
    _outputFile = new QLineEdit( ui->settingsFrame );
    _outputFile->setPlaceholderText( "result" );
    ui->settingsFrame->layout()->addWidget( _outputFile );

    QHBoxLayout* logHB = new QHBoxLayout( ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( new QLabel( "Log directory:", ui->settingsFrame ) );
    _logDir = new QLineEdit( ui->settingsFrame );
    _logDir->setPlaceholderText( "/path/to/dir" );
    logHB->addWidget( _logDir );
    QPushButton* logDirButton = new QPushButton( "Browse", ui->settingsFrame );
    logHB->addWidget( logDirButton );
    ui->settingsFrame->layout()->addItem( logHB );
    connect( logDirButton, SIGNAL( clicked() ), this, SLOT( selectLogDir() ) );

    QHBoxLayout* meshHB = new QHBoxLayout( ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( new QLabel( "Mesh file:", ui->settingsFrame ) );
    _meshFile = new QLineEdit( ui->settingsFrame );
    _meshFile->setPlaceholderText( "/path/to/mesh.su2" );
    meshHB->addWidget( _meshFile );
    QPushButton* meshFileButton = new QPushButton( "Browse", ui->settingsFrame );
    meshHB->addWidget( meshFileButton );
    ui->settingsFrame->layout()->addItem( meshHB );
    connect( meshFileButton, SIGNAL( clicked() ), this, SLOT( selectMeshFile() ) );

    ui->settingsFrame->layout()->addItem( new QSpacerItem( 1, 10, QSizePolicy::Fixed, QSizePolicy::Fixed ) );
    QFrame* hline = new QFrame( ui->settingsFrame );
    hline->setFrameShape( QFrame::HLine );
    hline->setFrameShadow( QFrame::Sunken );
    ui->settingsFrame->layout()->addWidget( hline );
    ui->settingsFrame->layout()->addItem( new QSpacerItem( 1, 10, QSizePolicy::Fixed, QSizePolicy::Fixed ) );

    //--------------------------------------------------

    ui->settingsFrame->layout()->addWidget( new QLabel( "Problem:", ui->settingsFrame ) );
    _problem = new QComboBox( ui->settingsFrame );
    for( auto& i : Problem_Map ) {
        _problem->addItem( QString::fromStdString( i.first ) );
    }
    ui->settingsFrame->layout()->addWidget( _problem );

    ui->settingsFrame->layout()->addWidget( new QLabel( "Solver:", ui->settingsFrame ) );
    _solver = new QComboBox( ui->settingsFrame );
    for( auto& i : Solver_Map ) {
        _solver->addItem( QString::fromStdString( i.first ) );
    }
    ui->settingsFrame->layout()->addWidget( _solver );

    ui->settingsFrame->layout()->addWidget( new QLabel( "CFL number:", ui->settingsFrame ) );
    _cfl = new QLineEdit( ui->settingsFrame );
    _cfl->setValidator( new QDoubleValidator( 0, 10, 10, this ) );
    _cfl->setPlaceholderText( "e.g. 0.9" );
    ui->settingsFrame->layout()->addWidget( _cfl );

    ui->settingsFrame->layout()->addWidget( new QLabel( "Simulation end time:", ui->settingsFrame ) );
    _tEnd = new QLineEdit( ui->settingsFrame );
    _tEnd->setValidator( new QDoubleValidator( 0, 10, 10, this ) );
    _tEnd->setPlaceholderText( "e.g. 1.0" );
    ui->settingsFrame->layout()->addWidget( _tEnd );

    QHBoxLayout* bcHB = new QHBoxLayout( ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( new QLabel( "Boundary Condition:", ui->settingsFrame ) );
    _bc = new QComboBox( ui->settingsFrame );
    _bc->addItem( QString::fromStdString( "BC_DIRICHLET" ) );
    _bc->addItem( QString::fromStdString( "BC_NEUMANN" ) );
    _bcNames = new QLineEdit( ui->settingsFrame );
    _bcNames->setPlaceholderText( "SU2 boundary names" );
    bcHB->addWidget( _bc );
    bcHB->addWidget( _bcNames );
    ui->settingsFrame->layout()->addItem( bcHB );

    ui->settingsFrame->layout()->addItem( new QSpacerItem( 1, 10, QSizePolicy::Fixed, QSizePolicy::Fixed ) );
    QFrame* hline2 = new QFrame( ui->settingsFrame );
    hline2->setFrameShape( QFrame::HLine );
    hline2->setFrameShadow( QFrame::Sunken );
    ui->settingsFrame->layout()->addWidget( hline2 );
    ui->settingsFrame->layout()->addItem( new QSpacerItem( 1, 10, QSizePolicy::Fixed, QSizePolicy::Fixed ) );

    //--------------------------------------------------

    ui->settingsFrame->layout()->addWidget( new QLabel( "Quadrature type:", ui->settingsFrame ) );
    _quadType = new QComboBox( ui->settingsFrame );
    for( auto& i : Quadrature_Map ) {
        _quadType->addItem( QString::fromStdString( i.first ) );
    }
    ui->settingsFrame->layout()->addWidget( _quadType );

    ui->settingsFrame->layout()->addWidget( new QLabel( "Quadrature order:", ui->settingsFrame ) );
    _quadOrder = new QLineEdit( ui->settingsFrame );
    _quadOrder->setValidator( new QIntValidator( 1, 1e6, this ) );
    _quadOrder->setPlaceholderText( "e.g. 8" );
    ui->settingsFrame->layout()->addWidget( _quadOrder );

    //--------------------------------------------------

    ui->settingsFrame->layout()->addItem( new QSpacerItem( 1, 1, QSizePolicy::Fixed, QSizePolicy::Expanding ) );

    QPushButton* load = new QPushButton( "Load", ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( load );
    load->setEnabled( false );    // todo: add functionality
    load->setStyleSheet( "color: gray" );
    QPushButton* save = new QPushButton( "Save", ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( save );
    QPushButton* run = new QPushButton( "Run", ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( run );

    //--------------------------------------------------

    _renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    _plotWindow   = new QVTKOpenGLWidget();
    _plotWindow->setMinimumSize( QSize( 800, 600 ) );
    _plotWindow->SetRenderWindow( _renderWindow.Get() );
    ui->plotTab->layout()->addWidget( _plotWindow );
    _renderer = vtkSmartPointer<vtkRenderer>::New();
    _plotWindow->GetRenderWindow()->AddRenderer( _renderer.Get() );

    ui->logViewer->setStyleSheet( "background-color: black" );
    ui->logViewer->setFont( QFontDatabase::systemFont( QFontDatabase::FixedFont ) );

    ui->tabWidget->setCurrentIndex( 0 );

    setStatusBarText( "Idle" );
    new QShortcut( QKeySequence( Qt::CTRL + Qt::Key_Q ), this, SLOT( close() ) );

    _logfileWatcher = new QFileSystemWatcher( this );

    connect( _logfileWatcher, SIGNAL( fileChanged( const QString& ) ), this, SLOT( appendLog( const QString& ) ) );
    // connect( load, SIGNAL( clicked() ), this, SLOT( /*TODO*/ ) );    // todo: add functionality
    connect( save, SIGNAL( clicked() ), this, SLOT( saveInputFile() ) );
    connect( run, SIGNAL( clicked() ), this, SLOT( runSimulation() ) );

    // quick debug settings
    //_outputDir->setText( "/home/jannick/Projects/rtsn/code/result" );
    //_outputFile->setText( "test" );
    //_logDir->setText( "/home/jannick/Projects/rtsn/code/result/logs" );
    //_meshFile->setText( "/home/jannick/Projects/rtsn/code/tests/input/mesh_files/checkerboard.su2" );
    //_problem->setCurrentIndex( 1 );
    //_solver->setCurrentIndex( 8 );
    //_cfl->setText( "0.5" );
    //_tEnd->setText( "0.4" );
    //_bc->setCurrentIndex( 1 );
    //_bcNames->setText( "void" );
    //_quadType->setCurrentIndex( 3 );
    //_quadOrder->setText( "15" );
}

void MainWindow::setStatusBarText( const QString& text ) { statusBar()->showMessage( text ); }

void MainWindow::appendLog( const QString& logfile ) {
    ui->logViewer->clear();
    QFile file( logfile );
    if( file.open( QIODevice::ReadOnly ) ) {

        ui->logViewer->appendPlainText( file.readAll() );
    }
    QCoreApplication::processEvents();
    this->update();
}

void MainWindow::runSimulation() {
    setStatusBarText( "Simulating..." );
    std::string tmpFile = std::filesystem::temp_directory_path().append( "KiT-RT_gui.cfg" );
    writeInputFile( QString::fromStdString( tmpFile ) );
    ui->tabWidget->setCurrentIndex( 1 );
    QCoreApplication::processEvents();
    Config* config = new Config( tmpFile );
    _logfileWatcher->addPath( QString::fromStdString( config->GetLogDir() + config->GetLogFile() ) );
    PrintLogHeader( tmpFile );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->PrintVolumeOutput();
    ui->tabWidget->setCurrentIndex( 0 );
    QString vtkFile = QString::fromStdString( config->GetOutputFile() );
    if( !vtkFile.endsWith( ".vtk" ) ) vtkFile.append( ".vtk" );
    plotResult( vtkFile );
    delete solver;
    delete config;
    std::filesystem::remove( tmpFile );
    setStatusBarText( "Idle" );
}

void MainWindow::plotResult( const QString filepath ) {
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName( filepath.toStdString().c_str() );
    reader->ReadAllScalarsOn();
    reader->ReadAllVectorsOn();
    reader->Update();
    auto output = reader->GetOutput();

    auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    ctf->SetColorSpaceToRGB();
    ctf->AddRGBPoint( -1.0, 1.0, 0.0, 0.0 );
    for( float i = 0; i < viridis.rows(); ++i ) {
        ctf->AddRGBPoint( viridis( i, 0 ), viridis( i, 1 ), viridis( i, 2 ), viridis( i, 3 ) );
    }
    auto lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetNumberOfTableValues( 256 );
    lut->Build();
    for( unsigned i = 0; i < 256; ++i ) {
        auto rgb = ctf->GetColor( i / 255.0 );
        lut->SetTableValue( i, rgb[0], rgb[1], rgb[2] );
    }

    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData( output );
    mapper->ScalarVisibilityOn();
    mapper->SetColorModeToMapScalars();
    mapper->SetLookupTable( lut );
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray( 0 );
    mapper->SetScalarRange( output->GetCellData()->GetArray( 0 )->GetRange() );

    auto scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable( mapper->GetLookupTable() );
    scalarBar->SetTitle( output->GetCellData()->GetArrayName( 0 ) );
    scalarBar->SetOrientationToHorizontal();
    scalarBar->SetPosition( 0.1, -0.001 );
    scalarBar->SetLabelFormat( "%-#6.1e" );
    scalarBar->SetWidth( 0.8 );
    scalarBar->SetHeight( 0.1 );
    scalarBar->SetNumberOfLabels( 4 );
    scalarBar->SetMaximumNumberOfColors( 256 );
    scalarBar->SetTitleRatio( 0.6 );
    auto titleprop = scalarBar->GetTitleTextProperty();
    titleprop->ShadowOff();
    titleprop->BoldOff();
    titleprop->SetColor( 1, 1, 1 );
    auto labelprop = scalarBar->GetLabelTextProperty();
    labelprop->ShadowOff();
    labelprop->BoldOff();
    labelprop->SetColor( 1, 1, 1 );
    scalarBar->SetLabelTextProperty( labelprop );

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper( mapper );

    _renderer->AddActor( actor );
    _renderer->AddActor2D( scalarBar );
    _renderer->ResetCamera();
    _renderer->UseFXAAOn();

    _plotWindow->update();
}

void MainWindow::writeInputFile( const QString filename ) {
    std::ofstream file( filename.toStdString() );
    file << "% AUTOGENERATED BY KIT-RT_GUI\n" << std::endl;
    file << "OUTPUT_DIR = " << _outputDir->text().toStdString() << std::endl;
    file << "OUTPUT_FILE = " << _outputFile->text().toStdString() << std::endl;
    file << "LOG_DIR = " << _logDir->text().toStdString() << std::endl;
    file << "MESH_FILE = " << _meshFile->text().toStdString() << std::endl;
    file << "PROBLEM = " << _problem->currentText().toStdString() << std::endl;
    file << "CFL_NUMBER = " << _cfl->text().toStdString() << std::endl;
    file << "TIME_FINAL = " << _tEnd->text().toStdString() << std::endl;
    file << _bc->currentText().toStdString() << " = ( " << _bcNames->text().toStdString() << " )" << std::endl;
    file << "SOLVER = " << _solver->currentText().toStdString() << std::endl;
    file << "QUAD_TYPE = " << _quadType->currentText().toStdString() << std::endl;
    file << "QUAD_ORDER = " << _quadOrder->text().toStdString() << std::endl;
    file.close();
}

void MainWindow::selectOutputDir() {
    QString dirName = QFileDialog::getExistingDirectory( this,
                                                         tr( "Select Output Directory" ),
                                                         QString::fromStdString( std::filesystem::current_path().string() ),
                                                         QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks );
    _outputDir->setText( dirName );
}
void MainWindow::selectLogDir() {
    QString dirName = QFileDialog::getExistingDirectory( this,
                                                         tr( "Select Log Directory" ),
                                                         QString::fromStdString( std::filesystem::current_path().string() ),
                                                         QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks );
    _logDir->setText( dirName );
}

void MainWindow::selectMeshFile() {
    QString fileName = QFileDialog::getOpenFileName( this, tr( "Load Mesh File" ), "", tr( "SU2 Mesh File (*.su2)" ) );
    _meshFile->setText( fileName );
}

void MainWindow::saveInputFile() {
    QString fileName = QFileDialog::getSaveFileName( this, tr( "Save Config File" ), "", tr( "Config File (*.cfg)" ) );
    if( !fileName.endsWith( ".cfg" ) ) fileName.append( ".cfg" );
    writeInputFile( fileName );
}
