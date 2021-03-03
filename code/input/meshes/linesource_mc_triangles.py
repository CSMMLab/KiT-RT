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

char_length = 0.01
geom = pg.opencascade.Geometry()
domain = add_block(0, 0, 0.99, 0.99, char_length, geom)
source = add_block(0.49, 0.49, 0.01, 0.01, char_length, geom)
geom.boolean_fragments([domain, source],[])
geom.add_physical(domain.lines, label="void")

mesh_code = geom.get_code()
with open("linesource_mc2.geo","w") as mesh_file:
    mesh_file.write(mesh_code)
os.system('gmsh linesource_mc2.geo -2 -format su2 -save_all')
os.system('rm linesource_mc2.geo')
