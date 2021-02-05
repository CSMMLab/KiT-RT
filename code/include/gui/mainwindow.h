#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#define vtkRenderingCore_AUTOINIT 2( vtkRenderingOpenGL2, vtkInteractionStyle )

#include <QComboBox>
#include <QFileSystemWatcher>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QPushButton>
#include <QShortcut>
#include <QString>

#include <QVTKOpenGLWidget.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkColorTransferFunction.h>
#include <vtkDataSetMapper.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkLookupTable.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

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
    QFileSystemWatcher* _logfileWatcher;
    QVTKOpenGLWidget* _plotWindow;
    vtkSmartPointer<vtkRenderer> _renderer;
    vtkSmartPointer<vtkRenderWindow> _renderWindow;

    void initUI();
    void plotResult( const QString filepath );

  public:
    explicit MainWindow();
    virtual ~MainWindow();

  public slots:
    void setStatusBarText( const QString& text );
    void appendLog( const QString& logfile );
    void runSimulation();

  private:
    Ui::MainWindow* ui;
};

#endif    // MAINWINDOW_H
