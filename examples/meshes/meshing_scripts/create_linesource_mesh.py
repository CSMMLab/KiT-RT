import pygmsh as pg
import numpy as np
import itertools
import os
from optparse import OptionParser

    
def add_block(x0,y0,lengthX, lengthY,char_length,geom):
    coords = np.array([
        [x0, y0, 0.0],
        [x0+lengthX, y0, 0.0],
        [x0+lengthX, y0+lengthY, 0.0],
        [x0, y0+lengthY, 0.0]
    ])
    return geom.add_polygon(coords, char_length)


def main():
    print("---------- Start Creating the mesh ------------")
    print("Parsing options")
    # --- parse options ---
    parser = OptionParser()
    parser.add_option("-o", "--output_name", dest="output_name", default="linesource")
    parser.add_option("-c", "--char_length", dest="char_length", default= 0.015)
    (options, args) = parser.parse_args()

    options.output_name = str(options.output_name)
    options.char_length = float(options.char_length)
    char_length = options.char_length
    geom = pg.opencascade.Geometry()
    domain = add_block(-1, -1, 2,2, char_length, geom)
    #geom.add_raw_code('psource = newp;\nPoint(psource) = {0.0, 0.0, 0.0, '+str(char_length)+'};\nPoint{psource} In Surface{'+domain.id+'};')
    geom.add_physical(domain.lines, label="void")


    mesh_code = geom.get_code()
    with open( options.output_name +".geo","w") as mesh_file:
        mesh_file.write(mesh_code)
    os.system('gmsh ' + options.output_name +'.geo -2 -format su2 -save_all')
    #os.system('rm ' + options.output_name +'.geo')
    print("---------- Successfully created the mesh ------------")

if __name__ == '__main__':
    main()
