import pygmsh as pg
import numpy as np
import itertools
import os
from optparse import OptionParser

    
def add_block(x0,y0,length,char_length,geom):
    coords = np.array([
        [x0, y0, 0.0],
        [x0+length, y0, 0.0],
        [x0+length, y0+length, 0.0],
        [x0, y0+length, 0.0]
    ])
    return geom.add_polygon(coords, char_length)



def main():
    print("---------- Start Creating the mesh ------------")
    print("Parsing options")
    # --- parse options ---
    parser = OptionParser()
    parser.add_option("-o", "--output_name", dest="output_name", default="lattice_unstructured")
    parser.add_option("-c", "--char_length", dest="char_length", default= 0.075)
    (options, args) = parser.parse_args()

    options.output_name = str(options.output_name)
    options.char_length = float(options.char_length)
    char_length = options.char_length
    geom = pg.opencascade.Geometry()
    domain = add_block(-3.5, -3.5, 7, char_length, geom)
    xpos = ypos = [1-3.5, 2-3.5, 3-3.5, 4-3.5, 5-3.5]
    pattern = list(itertools.product(xpos, ypos))[::2]
    pattern.pop(7)
    boxes = [domain]
    for pos in pattern:
        boxes.append(add_block(pos[0], pos[1], 1, char_length, geom))
    geom.boolean_fragments(boxes,[])
    geom.add_physical(domain.lines, label="void")

    mesh_code = geom.get_code()
    with open( options.output_name +".geo","w") as mesh_file:
        mesh_file.write(mesh_code)
    os.system('gmsh ' + options.output_name +'.geo -2 -format su2 -save_all')
    os.system('rm ' + options.output_name +'.geo')
    print("---------- Successfully created the mesh ------------")

if __name__ == '__main__':
    main()
