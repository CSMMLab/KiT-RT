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
    parser.add_option("-o", "--output_name", dest="output_name", default="struct_mesh_symmetric_hohlraum.su2")
    parser.add_option("-c", "--char_length", dest="char_length", default=0.01)
    parser.add_option("-s", "--start_pt", dest="start_pt", nargs=2, default=(-0.65,-0.65))
    parser.add_option("-l", "--length", dest="length", nargs=2, default=(1.3,1.3))
    parser.add_option("-b", "--boundary", dest="b_type", default="dirichletNeumann")
    (options, args) = parser.parse_args()

    options.output_name = str(options.output_name)
    options.char_length = float(options.char_length)
    options.b_type = str(options.b_type)
    char_length = options.char_length
    lengthX = float(options.length[0])
    lengthY = float(options.length[1])
    x0 = float(options.start_pt[0])
    y0 = float(options.start_pt[1])

    geom = pg.opencascade.Geometry()
    domain = add_block(x0,y0,lengthX, lengthY,char_length,geom)

    if options.b_type == "void":
        geom.add_physical(domain.lines, label="void")

    elif options.b_type == "voidAdiabatic":
        adiabatic = list()
        void = list()
        adiabatic.append(domain.lines[0])
        adiabatic.append(domain.lines[2])

        void.append(domain.lines[1])
        void.append(domain.lines[3])

        geom.add_physical(adiabatic, label="adiabatic")
        geom.add_physical(void, label="void")

    elif options.b_type == "dirichletNeumann":
        dirichlet_wall_up = list()
        dirichlet_wall_low = list()
        neumann_wall_left = list()
        neumann_wall_right = list()
        dirichlet_wall_low.append(domain.lines[0])
        dirichlet_wall_up.append(domain.lines[2])

        neumann_wall_left.append(domain.lines[1])
        neumann_wall_right.append(domain.lines[3])

        geom.add_physical(dirichlet_wall_up, label="dirichlet_wall_low")
        geom.add_physical(dirichlet_wall_low, label="dirichlet_wall_up")
        geom.add_physical(neumann_wall_left, label="neumann_wall_left")
        geom.add_physical(neumann_wall_right, label="neumann_wall_right")
        print("| Neumann marker has tag neumann_wall_left and neumann_wall_right")
        print("| Dirichlet marker has tag dirichlet_wall_up and dirichlet_wall_low")

    else:
        print("Boundary condition not yet implemented")

    geom.add_raw_code("Transfinite Surface{:};")
    geom.add_raw_code('Mesh.Algorithm= 3;')
    geom.add_raw_code("Recombine Surface {:} = 0;")

    mesh_code = geom.get_code()
    with open(options.output_name + ".geo", "w") as mesh_file:
        mesh_file.write(mesh_code)
    os.system('gmsh ' + options.output_name + '.geo -2 -format su2 -save_all')
    os.system('rm ' + options.output_name + '.geo')
    print("---------- Successfully created the mesh ------------")


if __name__ == '__main__':
    main()
