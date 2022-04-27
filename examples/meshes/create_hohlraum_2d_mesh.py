import pygmsh as pg
import numpy as np
import itertools
import os
from optparse import OptionParser


def add_block(x0, y0, x_len, y_len, char_length, geom):
    coords = np.array([
        [x0, y0, 0.0],
        [x0 + x_len, y0, 0.0],
        [x0 + x_len, y0 + y_len, 0.0],
        [x0, y0 + y_len, 0.0]
    ])
    return geom.add_polygon(coords, char_length)


def main():
    print("---------- Start Creating the mesh ------------")
    print("Parsing options")
    # --- parse options ---
    parser = OptionParser()
    parser.add_option("-o", "--output_name", dest="output_name", default="hohlraum")
    parser.add_option("-c", "--char_length", dest="char_length", default=0.008)
    (options, args) = parser.parse_args()

    options.output_name = str(options.output_name)
    options.char_length = float(options.char_length)
    char_length = options.char_length
    geom = pg.opencascade.Geometry()
    domain = add_block(0, 0, 1.3, 1.3, char_length, geom)
    xpos = [0, 0.45, 0.5, 0.5, 0.5, 0.0, 0.0, 1.25]
    ypos = [0.25, 0.25, 0.25, 1.0, 0.3, 0.0, 1.25, 0.05]
    x_len = [0.05, 0.05, 0.35, 0.35, 0.35, 1.3, 1.3, 0.05]
    y_len = [0.8, 0.8, 0.05, 0.05, 0.7, 0.05, 0.05, 1.2]

    pattern = list(itertools.product(xpos, ypos))[::2]
    pattern.pop(7)
    boxes = [domain]
    for i in range(len(xpos)):
        boxes.append(add_block(xpos[i], ypos[i], x_len[i], y_len[i], char_length, geom))
    geom.boolean_fragments(boxes, [])
    geom.add_physical(domain.lines, label="void")

    # geom.add_raw_code("Transfinite Surface{:};")
    # geom.add_raw_code('Mesh.Algorithm= 3;')
    # geom.add_raw_code("Recombine Surface {:} = 0;")

    mesh_code = geom.get_code()
    with open(options.output_name + ".geo", "w") as mesh_file:
        mesh_file.write(mesh_code)
    os.system('gmsh ' + options.output_name + '.geo -2 -format su2 -save_all')
    os.system('rm ' + options.output_name + '.geo')
    print("---------- Successfully created the mesh ------------")


if __name__ == '__main__':
    main()
