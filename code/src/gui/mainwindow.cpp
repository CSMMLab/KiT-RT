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
    QLineEdit* outputDir = new QLineEdit( ui->settingsFrame );
    outputDir->setPlaceholderText( "/path/to/dir" );
    outputHB->addWidget( outputDir );
    QPushButton* outputDirButton = new QPushButton( "Browse", ui->settingsFrame );
    outputHB->addWidget( outputDirButton );
    ui->settingsFrame->layout()->addItem( outputHB );

    ui->settingsFrame->layout()->addWidget( new QLabel( "Output filename:", ui->settingsFrame ) );
    QLineEdit* outputFile = new QLineEdit( ui->settingsFrame );
    outputFile->setPlaceholderText( "result" );
    ui->settingsFrame->layout()->addWidget( outputFile );

    QHBoxLayout* logHB = new QHBoxLayout( ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( new QLabel( "Log directory:", ui->settingsFrame ) );
    QLineEdit* logDir = new QLineEdit( ui->settingsFrame );
    logDir->setPlaceholderText( "/path/to/dir" );
    logHB->addWidget( logDir );
    QPushButton* logDirButton = new QPushButton( "Browse", ui->settingsFrame );
    logHB->addWidget( logDirButton );
    ui->settingsFrame->layout()->addItem( logHB );

    QHBoxLayout* meshHB = new QHBoxLayout( ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( new QLabel( "Mesh file:", ui->settingsFrame ) );
    QLineEdit* meshFile = new QLineEdit( ui->settingsFrame );
    meshFile->setPlaceholderText( "/path/to/mesh.su2" );
    meshHB->addWidget( meshFile );
    QPushButton* meshFileButton = new QPushButton( "Browse", ui->settingsFrame );
    meshHB->addWidget( meshFileButton );
    ui->settingsFrame->layout()->addItem( meshHB );

    ui->settingsFrame->layout()->addItem( new QSpacerItem( 1, 10, QSizePolicy::Fixed, QSizePolicy::Fixed ) );
    QFrame* hline = new QFrame( ui->settingsFrame );
    hline->setFrameShape( QFrame::HLine );
    hline->setFrameShadow( QFrame::Sunken );
    ui->settingsFrame->layout()->addWidget( hline );
    ui->settingsFrame->layout()->addItem( new QSpacerItem( 1, 10, QSizePolicy::Fixed, QSizePolicy::Fixed ) );

    //--------------------------------------------------

    ui->settingsFrame->layout()->addWidget( new QLabel( "Problem:", ui->settingsFrame ) );
    QComboBox* problem = new QComboBox( ui->settingsFrame );
    for( auto& i : Problem_Map ) {
        problem->addItem( QString::fromStdString( i.first ) );
    }
    ui->settingsFrame->layout()->addWidget( problem );

    ui->settingsFrame->layout()->addWidget( new QLabel( "Solver:", ui->settingsFrame ) );
    QComboBox* solver = new QComboBox( ui->settingsFrame );
    for( auto& i : Solver_Map ) {
        solver->addItem( QString::fromStdString( i.first ) );
    }
    ui->settingsFrame->layout()->addWidget( solver );

    ui->settingsFrame->layout()->addWidget( new QLabel( "CFL number:", ui->settingsFrame ) );
    QLineEdit* cfl = new QLineEdit( ui->settingsFrame );
    cfl->setValidator( new QDoubleValidator( 0, 10, 10, this ) );
    cfl->setPlaceholderText( "0.9" );
    ui->settingsFrame->layout()->addWidget( cfl );

    ui->settingsFrame->layout()->addWidget( new QLabel( "Simulation end time:", ui->settingsFrame ) );
    QLineEdit* tEnd = new QLineEdit( ui->settingsFrame );
    tEnd->setValidator( new QDoubleValidator( 0, 10, 10, this ) );
    tEnd->setPlaceholderText( "1.0" );
    ui->settingsFrame->layout()->addWidget( tEnd );

    ui->settingsFrame->layout()->addItem( new QSpacerItem( 1, 10, QSizePolicy::Fixed, QSizePolicy::Fixed ) );
    QFrame* hline2 = new QFrame( ui->settingsFrame );
    hline2->setFrameShape( QFrame::HLine );
    hline2->setFrameShadow( QFrame::Sunken );
    ui->settingsFrame->layout()->addWidget( hline2 );
    ui->settingsFrame->layout()->addItem( new QSpacerItem( 1, 10, QSizePolicy::Fixed, QSizePolicy::Fixed ) );

    //--------------------------------------------------

    ui->settingsFrame->layout()->addWidget( new QLabel( "Quadrature type:", ui->settingsFrame ) );
    QComboBox* quadType = new QComboBox( ui->settingsFrame );
    for( auto& i : Quadrature_Map ) {
        quadType->addItem( QString::fromStdString( i.first ) );
    }
    ui->settingsFrame->layout()->addWidget( quadType );

    ui->settingsFrame->layout()->addWidget( new QLabel( "Quadrature order:", ui->settingsFrame ) );
    QLineEdit* quadOrder = new QLineEdit( ui->settingsFrame );
    quadOrder->setValidator( new QIntValidator( 1, 1e6, this ) );
    quadOrder->setPlaceholderText( "8" );
    ui->settingsFrame->layout()->addWidget( quadOrder );

    //--------------------------------------------------

    ui->settingsFrame->layout()->addItem( new QSpacerItem( 1, 1, QSizePolicy::Fixed, QSizePolicy::Expanding ) );

    QPushButton* load = new QPushButton( "Load", ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( load );
    QPushButton* save = new QPushButton( "Save", ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( save );
    QPushButton* run = new QPushButton( "Run", ui->settingsFrame );
    ui->settingsFrame->layout()->addWidget( run );

    //--------------------------------------------------

    _renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    _plotWindow   = new QVTKOpenGLWidget();
    _plotWindow->SetRenderWindow( _renderWindow.Get() );
    ui->plotTab->layout()->addWidget( _plotWindow );
    _renderer = vtkSmartPointer<vtkRenderer>::New();
    _plotWindow->GetRenderWindow()->AddRenderer( _renderer.Get() );

    ui->tabWidget->setCurrentIndex( 0 );

    setStatusBarText( "Idle" );
    new QShortcut( QKeySequence( Qt::CTRL + Qt::Key_Q ), this, SLOT( close() ) );

    _logfileWatcher = new QFileSystemWatcher( this );

    connect( _logfileWatcher, SIGNAL( fileChanged( const QString& ) ), ui->logViewer, SLOT( appendLog( const QString& ) ) );
    connect( run, SIGNAL( clicked() ), this, SLOT( runSimulation() ) );
}

void MainWindow::setStatusBarText( const QString& text ) { statusBar()->showMessage( text ); }

void MainWindow::appendLog( const QString& logfile ) {
    QFile file( logfile );
    if( file.open( QIODevice::ReadOnly ) ) {
        qint64 num      = 10;
        qint64 fileSize = file.size();
        file.seek( fileSize - num );
        ui->logViewer->appendPlainText( file.read( num ) );
    }
}

void MainWindow::runSimulation() {
    ui->tabWidget->setCurrentIndex( 1 );
    std::string filename = "Config needs a new constructor";
    Config* config       = new Config( filename );
    _logfileWatcher->addPath( QString::fromStdString( config->GetLogDir() + config->GetLogFile() ) );
    PrintLogHeader( filename );
    Solver* solver = Solver::Create( config );
    solver->Solve();
    solver->PrintVolumeOutput();
    ui->tabWidget->setCurrentIndex( 0 );
    plotResult( QString::fromStdString( config->GetOutputDir() + config->GetOutputFile() ) );
}

void MainWindow::plotResult( const QString filepath ) {
    vtkNew<vtkUnstructuredGridReader> reader;
    reader->SetFileName( filepath.toStdString().c_str() );
    reader->ReadAllScalarsOn();
    reader->ReadAllVectorsOn();
    reader->Update();
    auto output = reader->GetOutput();

    vtkNew<vtkColorTransferFunction> ctf;
    ctf->SetColorSpaceToRGB();
    Matrix cmap( 3, 3 );
    ctf->AddRGBPoint( -1, 1, 0, 0 );
    for( float i = 0; i < cmap.columns(); ++i ) {
        ctf->AddRGBPoint( i / 256.0, cmap( i, 0 ), cmap( i, 1 ), cmap( i, 2 ) );
    }
    vtkNew<vtkLookupTable> lut;
    lut->SetNumberOfTableValues( cmap.columns() );
    lut->Build();
    for( float i = 0; i < cmap.columns(); ++i ) {
        auto rgb = ctf->GetColor( i / 256.0 );
        lut->SetTableValue( i, rgb );
    }

    vtkNew<vtkCamera> camera;

    double camPos[3] = { 0.0, 0.0, 0.0 };
    double camFP[3]  = { 0.0, 0.0, 0.0 };
    camera->SetPosition( camPos[0], camPos[1], camPos[2] );
    camera->SetFocalPoint( camFP[0], camFP[1], camFP[2] );

    vtkNew<vtkDataSetMapper> mapper;
    mapper->SetInputData( output );
    mapper->ScalarVisibilityOn();
    mapper->SetColorModeToMapScalars();
    mapper->SetLookupTable( lut );
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray( 0 );
    mapper->SetScalarRange( output->GetCellData()->GetArray( 0 )->GetRange() );

    vtkNew<vtkScalarBarActor> scalarBar;
    scalarBar->SetLookupTable( mapper->GetLookupTable() );
    scalarBar->SetTitle( "foo" );
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
    titleprop->SetColor( 0, 0, 0 );
    auto labelprop = scalarBar->GetLabelTextProperty();
    labelprop->ShadowOff();
    labelprop->BoldOff();
    labelprop->SetColor( 0, 0, 0 );
    scalarBar->SetLabelTextProperty( labelprop );

    vtkNew<vtkActor> actor;
    actor->SetMapper( mapper );

    _renderer->AddActor( actor );
    _renderer->AddActor2D( scalarBar );
    _renderer->UseFXAAOn();
    _renderer->SetBackground( 1, 1, 1 );
    _renderer->SetActiveCamera( camera );

    _renderWindow->Render();
}
