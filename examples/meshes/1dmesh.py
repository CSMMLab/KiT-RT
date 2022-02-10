import pygmsh as pg
import numpy as np
import itertools
import os
from optparse import OptionParser


def add_block(x0,lengthX,char_length,geom):
    # coords = np.array([
    #     [x0, 0.0, 0.0],
    #     [x0+lengthX, 0.0, 0.0],
    #     [x0+lengthX, 0.5*char_length, 0.0],
    #     [x0, 0.5*char_length, 0.0]
    # ])
    coords = np.array([
        [x0, 0.0, 0.0],
        [x0+lengthX, 0.0, 0.0],
        [x0+lengthX, 0.3, 0.0],
        [x0, 0.3, 0.0]
    ])
    print(len(coords))
    print(len([char_length, 1]))
    return geom.add_polygon(coords, char_length)


def main():
    print("---------- Start Creating the mesh ------------")
    print("Parsing options")
    # --- parse options ---
    parser = OptionParser()
    parser.add_option("-o", "--output_name", dest="output_name", default="1dmesh")
    parser.add_option("-c", "--char_length", dest="char_length", default= 0.01)
    parser.add_option("-s", "--start_pt", dest="start_pt", default= 0)
    parser.add_option("-l", "--x_length", dest="x_length", default= 1)
    parser.add_option("-b", "--boundary", dest="b_type", default= "dirichletNeumann")
    (options, args) = parser.parse_args()

    options.output_name = str(options.output_name)
    options.char_length = float(options.char_length)
    options.b_type = str(options.b_type)
    options.x_length = float(options.x_length)
    char_length = options.char_length
    lengthX = options.x_length
    start_pt = float(options.start_pt)

    geom = pg.opencascade.Geometry()
    domain = add_block(start_pt, lengthX, char_length, geom)

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
        dirichlet = list()
        wall_low = list()
        wall_up = list()
        wall_low.append(domain.lines[0])
        wall_up.append(domain.lines[2])

        dirichlet.append(domain.lines[1])
        dirichlet.append(domain.lines[3])

        geom.add_physical(dirichlet, label="dirichlet")
        geom.add_physical(wall_low, label="wall_low")
        geom.add_physical(wall_up, label="wall_up")

    else:
        print("Boundary condition not yet implemented")

    geom.add_raw_code("Transfinite Curve {2,4} = 1 Using Progression 1;")
    geom.add_raw_code("Transfinite Surface{:};")
    geom.add_raw_code('Mesh.Algorithm= 3;')
    geom.add_raw_code("Recombine Surface {:} = 0;")
    

    mesh_code = geom.get_code()
    with open( options.output_name +".geo","w") as mesh_file:
        mesh_file.write(mesh_code)
    os.system('gmsh ' + options.output_name +'.geo -2 -format su2 -save_all')
    os.system('rm ' + options.output_name +'.geo')
    print("---------- Successfully created the mesh ------------")


if __name__ == '__main__':
    main()
