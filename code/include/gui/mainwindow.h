#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#define vtkRenderingCore_AUTOINIT 3( vtkRenderingOpenGL2, vtkInteractionStyle, vtkRenderingFreeType )

#include <QComboBox>
#include <QFileDialog>
#include <QFileSystemWatcher>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QPushButton>
#include <QShortcut>
#include <QString>
#include <QSurfaceFormat>

#include <QVTKOpenGLWidget.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkColorTransferFunction.h>
#include <vtkDataSetMapper.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkLookupTable.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

#include "cmap.h"
#include "common/config.h"
#include "common/globalconstants.h"
#include "common/io.h"
#include "solvers/solverbase.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
  private:
    Q_OBJECT

    QLineEdit* _outputDir;
    QLineEdit* _outputFile;
    QLineEdit* _logDir;
    QLineEdit* _meshFile;
    QLineEdit* _cfl;
    QLineEdit* _tEnd;
    QLineEdit* _bcNames;
    QLineEdit* _quadOrder;

    QComboBox* _problem;
    QComboBox* _solver;
    QComboBox* _bc;
    QComboBox* _quadType;

    QFileSystemWatcher* _logfileWatcher;
    QVTKOpenGLWidget* _plotWindow;
    vtkSmartPointer<vtkRenderer> _renderer;
    vtkSmartPointer<vtkRenderWindow> _renderWindow;

    void initUI();
    void plotResult( const QString filepath );
    void writeInputFile( const QString filename );

  public:
    explicit MainWindow();
    virtual ~MainWindow();

  public slots:
    void setStatusBarText( const QString& text );
    void appendLog( const QString& logfile );
    void runSimulation();

    void selectOutputDir();
    void selectLogDir();
    void selectMeshFile();

    void saveInputFile();

  private:
    Ui::MainWindow* ui;
};

#endif    // MAINWINDOW_H
