import matplotlib.pyplot as plt
from matplotlib import cm
import argparse
import numpy as np
import datetime as dt
import os
import pandas as pd
import re
import errno
import vtk
import itertools
import seaborn as sns

#################### SETUP ####################

def getProblemVTKSettings(problem):
    if problem == 'CHECKERBOARD':
        return {
            'logScale': True,
            'cameraPosition': [3.5, 3, -17.5],
            'cameraFocalPoint': [3.5, 3, 0]
        }
    else:
        return{
            'logScale': False,
            'cameraPosition': [0,0,0],
            'cameraFocalPoint': [0,0,0]        
        }

def getLabel(s):
    labelDict = {
        'Iter': 'Iterations',
        'Runtime': 'Runtime [s]',
        'Mass': 'Mass',
        'RMS_flux': 'L$_2$ flux difference',        
    }
    if s in labelDict:
        return labelDict[s]
    else:
        return s        

def latexify(s):
    if any(latexExpression in s for latexExpression in ['ρ','µ','θ']):
        s = ''.join('${0}$'.format(s))
        s.replace('ρ', '\\rho ')
        s.replace('µ', '\\mu ')
        s.replace('θ', '\\theta ')    
    return s

#################### CODE ####################    

def timeToFloat(time):
    t = dt.datetime.strptime(time.strip(), "%Y-%m-%d %H:%M:%S.%f")
    epoch = dt.datetime.utcfromtimestamp(0)
    return (t-epoch).total_seconds()    

class LogFile:
    def __init__(self, filepath):
        self.filepath = os.path.normpath(filepath)
        if os.path.isfile(self.filepath):
            with open(self.filepath, 'r') as content:
                print("Parsing:\t" + self.filepath)
                self.raw = content.readlines()  
        else:
            raise FileNotFoundError(self.filepath)
        self.outputDir = os.path.dirname(filepath)              
        self.logFile = os.path.basename(filepath)        
        self.header = []
        sectionCtr = 0 # 1 = banner, 2 = copyright, 3 = git/mpi info, 4 = config file, 5 = spacer, 6 = log
        for line in self.raw:
            if '----------' in line:
                sectionCtr += 1
                continue
            elif sectionCtr == 4:                
                self.header.append(line.split('|')[1].strip()) #strip time stamp and leading/trailing whitespaces
            elif sectionCtr > 4:
                break
        self.header = list(filter(None, self.header)) #strip empty lines                
        for line in self.header:        
            if 'MESH_FILE' in line: self.meshFile = line.split('=')[1].strip()
            if 'PROBLEM' in line: self.problem = line.split('=')[1].strip()
            if 'SOLVER' in line: self.solver = line.split('=')[1].strip() 
            if 'CFL_NUMBER' in line: self.cfl = float(line.split('=')[1].strip())
            if 'TIME_FINAL' in line: self.tEnd = float(line.split('=')[1].strip())
            if 'QUAD_TYPE' in line: self.quadType = line.split('=')[1].strip() 
            if 'QUAD_ORDER' in line: self.quadOrder = int(line.split('=')[1].strip())        
            if 'OUTPUT_FILE' in line: self.vtkFile = line.split('=')[1].strip()  
        if os.path.isfile(self.filepath + '_csv'):
            df = pd.read_csv(self.filepath + '_csv')       
            df.rename(columns={ df.columns[0]: "Time" }, inplace = True)
            t0 = timeToFloat(df['Time'][0])
            df['Runtime'] = np.zeros(len(df['Time']))
            for index, _ in df.iterrows():
                df.loc[index,'Runtime'] = timeToFloat(df['Time'][index]) - t0
            df = df.drop(columns=['Time', 'VTK_out', 'CSV_out'])
            df.columns = df.columns.str.strip()
            self.data = df
        if not self.vtkFile.endswith('.vtk'):
            self.vtkFile += '.vtk'   

def createVTKPlots(dir, logFiles, disableRescaling=False, showMesh=False, showIsolines=False):
    outputdir = 'plots/vtk'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    if not all(l.problem == logFiles[0].problem for l in logFiles):
        raise(ValueError)    
    problem = logFiles[0].problem
    setup = getProblemVTKSettings(problem)
    minVal = {}
    maxVal = {}
    fieldNames = []
    for log in logFiles:
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(log.vtkFile)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        output = reader.GetOutput()        
        for i in range(output.GetCellData().GetNumberOfArrays()):
            f = output.GetCellData().GetArrayName(i)
            if f not in fieldNames:
                fieldNames.append(f)            
            valRange = output.GetCellData().GetArray(i).GetRange()                      
            minVal[f] = valRange[0] #np.maximum(minVal[f] if f in minVal else -np.inf, valRange[0])
            maxVal[f] = valRange[1] #np.minimum(maxVal[f] if f in maxVal else  np.inf, valRange[1])

    ctf = vtk.vtkColorTransferFunction()
    ctf.SetColorSpaceToRGB()
    cmap = cm.get_cmap('viridis', 256)
    ctf.AddRGBPoint(-1, 1, 0, 0)
    for i,rgb in enumerate(cmap.colors):        
        ctf.AddRGBPoint(float(i/256),rgb[0],rgb[1],rgb[2])        

    for log in logFiles:
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(log.vtkFile)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        output = reader.GetOutput()                


        if setup['logScale']:        
            lut = vtk.vtkLogLookupTable()
        else: 
            lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(256)
        lut.Build()
        for i in range(0,256):
            rgb = list(ctf.GetColor(float(i)/256))+[1]
            lut.SetTableValue(i,rgb)

        camera = vtk.vtkCamera()

        camPos =  setup['cameraPosition']
        camFP =  setup['cameraFocalPoint']
        camera.SetPosition(camPos[0],camPos[1],camPos[2]);
        camera.SetFocalPoint(camFP[0],camFP[1],camFP[2]);
        
        for f in fieldNames:
            mapper = vtk.vtkDataSetMapper()
            mapper.SetInputData(output)
            mapper.ScalarVisibilityOn()
            mapper.SetColorModeToMapScalars()
            mapper.SetLookupTable(lut)
            mapper.SetScalarModeToUsePointFieldData()
            mapper.SelectColorArray(f)
            if disableRescaling:
                mapper.SetScalarRange(output.GetCellData().GetArray(f).GetRange())                
            else:
                mapper.SetScalarRange(minVal[f], maxVal[f])

            scalarBar = vtk.vtkScalarBarActor()
            scalarBar.SetLookupTable(mapper.GetLookupTable())
            scalarBar.SetTitle(latexify(f))
            scalarBar.SetOrientationToHorizontal()
            scalarBar.SetPosition(0.1,-0.001)
            scalarBar.SetLabelFormat('%-#6.1e')
            scalarBar.SetWidth(0.8)
            scalarBar.SetHeight(0.1)
            scalarBar.SetNumberOfLabels(4)
            scalarBar.SetMaximumNumberOfColors(256)
            scalarBar.SetTitleRatio(0.6)
            titleprop = scalarBar.GetTitleTextProperty()
            titleprop.ShadowOff()
            titleprop.BoldOff()
            titleprop.SetColor(0,0,0)
            labelprop = scalarBar.GetLabelTextProperty()
            labelprop.ShadowOff()
            labelprop.BoldOff()
            labelprop.SetColor(0,0,0)
            scalarBar.SetLabelTextProperty(labelprop)

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            if showMesh:
                actor.GetProperty().EdgeVisibilityOn()
                actor.GetProperty().SetLineWidth(1.0)

            if showIsolines:
                output.GetPointData().SetActiveScalars(f)
                contours = vtk.vtkContourFilter()
                contours.SetInputData(output)
                contours.GenerateValues(10, minVal[f], maxVal[f])
                
                contMapper = vtk.vtkPolyDataMapper()
                contMapper.SetInputConnection(contours.GetOutputPort())
                contMapper.ScalarVisibilityOff()
                
                contActor = vtk.vtkActor()
                contActor.SetMapper(contMapper) 
                #contActor.GetProperty().SetColor(1,1,1)                

            renderer = vtk.vtkRenderer()
            renderer.AddActor(actor)
            if showIsolines: renderer.AddActor(contActor)
            renderer.AddActor2D(scalarBar)
            renderer.UseFXAAOn()
            renderer.SetBackground(1, 1, 1)
            renderer.SetActiveCamera(camera)

            render_window = vtk.vtkRenderWindow()
            render_window.SetOffScreenRendering(True)
            render_window.AddRenderer(renderer)
            render_window.SetSize(800,1000)
            render_window.Render()

            windowToImageFilter = vtk.vtkWindowToImageFilter()
            windowToImageFilter.SetInput(render_window)
            windowToImageFilter.ReadFrontBufferOff()
            windowToImageFilter.Update()

            writer = vtk.vtkPNGWriter()
            writer.SetFileName(outputdir + '/' + log.logFile + "_" + f + ".png")
            writer.SetInputConnection(windowToImageFilter.GetOutputPort())
            writer.Write()
            print('Plot created:\t' + outputdir + '/' + log.logFile + "_" + f + ".png")   

def createConvergencePlots(logFiles):
    outputDir = 'plots/convergence'
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    columns = list(logFiles[0].data)
    if not all(list(l.data) == list(logFiles[0].data) for l in logFiles):
        raise(ValueError)    
    cpal = ["#886a4c","#7548cf","#5dbf48","#d150c6","#a1ae43","#8081d3","#d19b45","#8f4888","#73b990","#dc4d32","#6facc3","#cb4671","#4e6f3f","#cf9bad","#ab583e","#5b5f79"]
    xData = ['Iter', 'Runtime']
    yData = [y for y in columns if y not in xData]
    for x in xData:
        for y in yData:
            plt.cla()
            plt.clf()            
            palette = itertools.cycle(sns.color_palette(cpal))
            plt.xlabel(getLabel(x))
            plt.title(getLabel(y))                        
            for f in logFiles:                
                plt.semilogy(f.data.loc[:,x].values, f.data.loc[:,y].values, color=next(palette))           
            plotName = os.path.normpath(x + '_' + y +'.pdf')
            plt.savefig(os.path.join(outputDir, plotName))
            print('Plot created:\t' + os.path.join(outputDir, plotName))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    disableRescaling = showIsolines = showMesh = False
    parser.add_argument("--dir", "-d", type=str, required=False, default='logs')
    parser.add_argument("--disable-rescaling", "-r", action='store_true', dest='disableRescaling')
    parser.add_argument("--show-isolines", "-i", action='store_true', dest='showIsolines')
    parser.add_argument("--show-mesh", "-m", action='store_true', dest='showMesh')
    args = parser.parse_args()
    dir = os.path.normpath(args.dir)
    if not os.path.exists('plots'):
        os.makedirs('plots')
    files = os.listdir(dir)
    logFiles = []
    for f in files:
        if not f.endswith('_csv'):
            logFiles.append(LogFile(os.path.join(dir, f)))    
    createConvergencePlots(logFiles)
    createVTKPlots(os.getcwd(), logFiles, disableRescaling, showMesh, showIsolines)