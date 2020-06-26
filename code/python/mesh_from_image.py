import numpy as np
import pygmsh as pg
import os

import extract_grayscale_image as egs

def add_rectangle(x0,y0,width,height,char_length,geom):
    coords = np.array([
        [x0, y0, 0.0],
        [x0+width, y0, 0.0],
        [x0+width, y0+height, 0.0],
        [x0, y0+height, 0.0]
    ])
    return geom.add_polygon(coords, char_length)

def generate(image_name, mesh_name):
    if mesh_name.endswith('.su2'):
        mesh_name = os.path.splitext(mesh_name)[0]
    gsImage, dimensions = egs.extract(image_name)

    char_length = 1.0#np.min(1.0, np.min(dimensions)/10.0) #arbitrary
    geom = pg.opencascade.Geometry()
    domain = add_rectangle(0.0, 0.0, dimensions[0], dimensions[1], char_length, geom)
    geom.add_physical(domain.lines, label="void")

    mesh_code = geom.get_code()
    with open(mesh_name+".geo","w") as mesh_file:
        mesh_file.write(mesh_code,) 
    os.system('gmsh '+mesh_name+'.geo -2 -format su2 -save_all > /dev/null')
    os.system('rm '+mesh_name+'.geo > /dev/null')

    return gsImage
