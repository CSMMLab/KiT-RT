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
    parser.add_option("-c", "--char_length", dest="char_length", default=0.025)
    (options, args) = parser.parse_args()

    options.output_name = str(options.output_name)
    options.char_length = float(options.char_length)
    char_length = options.char_length
    geom = pg.opencascade.Geometry()
    domain = add_block(-0.05, -0.05, 1.4, 1.4, char_length, geom)

    boxes = [domain]

    # add left side absorption region
    boxes.append(add_block(0, 0.25, 0.05, 0.8, char_length, geom))

    # add interior polygon
    coords = np.array([
        [0.85, 0.25, 0.0],
        [0.45, 0.25, 0.0],
        [0.45, 1.05, 0.0],
        [0.85, 1.05, 0.0],
        [0.85, 1.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, 0.3, 0.0],
        [0.85, 0.3, 0.0],
    ])
    boxes.append(geom.add_polygon(coords, char_length))

    # add outer polygon
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.3, 0.0, 0.0],
        [1.3, 1.3, 0.0],
        [0.0, 1.3, 0.0],
        [0.0, 1.25, 0.0],
        [1.25, 1.25, 0.0],
        [1.25, 0.05, 0.0],
        [0.0, 0.05, 0.0],
    ])
    boxes.append(geom.add_polygon(coords, char_length))

    geom.boolean_fragments(boxes, [])
    geom.add_physical(domain.lines, label="void")

    # geom.add_raw_code("Transfinite Surface{:};")
    # geom.add_raw_code('Mesh.Algorithm= 3;')
    # geom.add_raw_code("Recombine Surface {:} = 0;")

    mesh_code = geom.get_code()
    with open(options.output_name + ".geo", "w") as mesh_file:
        mesh_file.write(mesh_code)
    os.system('gmsh ' + options.output_name + '.geo -2 -format su2 -save_all')
    # os.system('rm ' + options.output_name + '.geo')
    print("---------- Successfully created the mesh ------------")


if __name__ == '__main__':
    main()
