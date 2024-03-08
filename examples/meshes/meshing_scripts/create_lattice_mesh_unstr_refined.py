import pygmsh as pg
import numpy as np
import itertools
import os
from optparse import OptionParser
import subprocess


def add_block(x0, y0, length, char_length, geom):
    coords = np.array([
        [x0, y0, 0.0],
        [x0 + length, y0, 0.0],
        [x0 + length, y0 + length, 0.0],
        [x0, y0 + length, 0.0]
    ])
    return geom.add_polygon(coords, char_length)


def add_boundary(x0, y0, length, cl_max, cl_min, geom):
    coords = np.array([
        [x0, y0, 0.0],
        [x0 + length, y0, 0.0],
        [x0 + length, y0 + length, 0.0],
        [x0, y0 + length, 0.0]
    ])
    pts = []
    for coord in coords:
        pts.append(geom.add_point(coord, cl_max))

    l1 = geom.add_line(pts[0], pts[1])
    l2 = geom.add_line(pts[1], pts[2])
    l3 = geom.add_line(pts[2], pts[3])
    l4 = geom.add_line(pts[3], pts[0])

    bd1 = geom.add_boundary_layer(
        edges_list=[l1, l2],
        hfar=cl_max,
        hwall_n=cl_min,
        ratio=1.1,
        thickness=0.1,
        anisomax=100.0,
    )

    return bd1


def main():
    print("---------- Start Creating the mesh ------------")
    print("Parsing options")
    # --- parse options ---
    parser = OptionParser()
    parser.add_option("-o", "--output_name", dest="output_name", default="lattice_unstructured")
    parser.add_option("-c", "--cl_min", dest="cl_min", default=1e-3)
    parser.add_option("-d", "--cl_max", dest="cl_max", default=5 * 1e-1)
    (options, args) = parser.parse_args()

    options.output_name = str(options.output_name)
    options.cl_max = float(options.cl_max)
    options.cl_min = float(options.cl_min)
    cl_max = options.cl_max
    cl_min = options.cl_min

    geom = pg.opencascade.Geometry(characteristic_length_min=None, characteristic_length_max=None)
    domain = add_block(-3.5, -3.5, 7, cl_max, geom)

    xpos = ypos = [1 - 3.5, 2 - 3.5, 3 - 3.5, 4 - 3.5, 5 - 3.5]
    pattern = list(itertools.product(xpos, ypos))[::2]
    pattern.pop(7)
    
    boxes = [domain]
    boundary_layers = []
    thickness = 0.1
    # for pos in pattern:
    #    boxes.append(add_block(pos[0], pos[1], 1, cl_max, geom))
    # boundary_layers.append(add_boundary(pos[0], pos[1], 1, cl_max, cl_min, geom))
    # poly = boxes[5]
    pt1 = geom.add_point([0, 0, 0], lcar=cl_min)
    pt2 = geom.add_point([1, 1, 0], lcar=cl_min)
    pt3 = geom.add_point([0, 1, 0], lcar=cl_min)
    l0 = geom.add_line(pt1, pt2)

    # l1 = geom.add_line(pt2, pt3)
    # l2 = geom.add_line(pt3, pt1)
    # ll0 = geom.add_line_loop((l0, l1, l2))
    # rs0 = geom.add_surface(ll0)
    # tf_line = geom.set_transfinite_lines([l0], 100, progression=1.5)
    # geom.set_recombined_surfaces([rs0])

    # geom.boolean_fragments(boxes, [])

    geom.add_physical(domain.lines, label="void")

    mesh_code = geom.get_code()

    with open(options.output_name + ".geo", "w") as mesh_file:
        mesh_file.write(mesh_code)

    os.system('gmsh ' + options.output_name + '.geo -2 -format su2 -save_all')
    os.system('gmsh ' + options.output_name + '.geo -2 -format vtk -save_all')

    # os.system('rm ' + options.output_name + '.geo')
    # add mesh refinements
    print("---------- Successfully created the mesh ------------")


if __name__ == '__main__':
    main()
