import pygmsh as pg
import numpy as np
import itertools
import os

def add_block(x0,y0,length,char_length,geom):
    coords = np.array([
        [x0, y0, 0.0],
        [x0+length, y0, 0.0],
        [x0+length, y0+length, 0.0],
        [x0, y0+length, 0.0]
    ])
    return geom.add_polygon(coords, char_length)

char_length = 0.2
geom = pg.opencascade.Geometry()
domain = add_block(-0.6, -0.6, 1.2, char_length, geom)
geom.add_raw_code('psource = newp;\nPoint(psource) = {0.0, 0.0, 0.0, '+str(char_length)+'};\nPoint{psource} In Surface{'+domain.id+'};')
geom.add_physical(domain.lines, label="void")

mesh_code = geom.get_code()
with open("linesource.geo","w") as mesh_file:
    mesh_file.write(mesh_code)
os.system('gmsh linesource.geo -2 -format su2 -save_all')
os.system('rm linesource.geo')
