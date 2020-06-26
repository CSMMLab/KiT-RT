import pygmsh as pg
import numpy as np
import itertools
import os

def add_block(x0,y0,lengthX, lengthY,char_length,geom):
    coords = np.array([
        [x0, y0, 0.0],
        [x0+lengthX, y0, 0.0],
        [x0+lengthX, y0+lengthY, 0.0],
        [x0, y0+lengthY, 0.0]
    ])
    return geom.add_polygon(coords, char_length)

char_length = 0.015
geom = pg.opencascade.Geometry()
domain = add_block(-10, -1, 20,2, char_length, geom)
#geom.add_raw_code('psource = newp;\nPoint(psource) = {0.0, 0.0, 0.0, '+str(char_length)+'};\nPoint{psource} In Surface{'+domain.id+'};')

geom.add_physical(domain.lines[0], label="adiabatic")
geom.add_physical(domain.lines[1], label="void")
geom.add_physical(domain.lines[2], label="adiabatic")
geom.add_physical(domain.lines[3], label="void")

mesh_code = geom.get_code()
with open("linesource_pseudo_1D.geo","w") as mesh_file:
    mesh_file.write(mesh_code)
os.system('gmsh linesource_pseudo_1D.geo -2 -format su2 -save_all')
os.system('rm linesource_pseudo_1D.geo')
