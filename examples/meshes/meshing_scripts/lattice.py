# ------------------------------------------------------------------------------
#
#  Gmsh Python tutorial 12
#
#  Cross-patch meshing with compounds
#
# ------------------------------------------------------------------------------

import gmsh
import sys
from optparse import OptionParser


def main(filename, lc, sigma_t):
    # Specify domain
    refinement_box_width = 0.02
    transition_thickness = 0.02

    x = [-3.5, 3.5]
    y = [-3.5, 3.5]

    x_inner = [-2.5, -2.5, -2.5, -1.5, -1.5, -0.5, -0.5, .5, .5, 1.5, 1.5, 1.5]
    y_inner = [1.5, -0.5, -2.5, 0.5, -1.5, -0.5, -2.5, 0.5, -1.5, 1.5, -0.5, -2.5, ]

    gmsh.initialize()  # initialize gmsh

    domain = gmsh.model.occ.addRectangle(x[0], y[0], 0, 7, 7)

    gmsh.model.occ.synchronize()

    geometries = [(2, domain)]
    for x_i, y_i in zip(x_inner, y_inner):
        box = gmsh.model.occ.addRectangle(x_i, y_i, 0, 1, 1)
        geometries.append((2, box))

    ov, ovv = gmsh.model.occ.fragment(geometries, [])  # stitch together geometries
    gmsh.model.occ.synchronize()

    # Generate refinement areas in mesh
    box_x = [(-2.5, -2.5), (-2.5, -2.5), (-2.5, -2.5), (-1.5, -1.5), (-.5, -.5), (0.5, 0.5), (1.5, 1.5), (2.5, 2.5),
             (2.5, 2.5), (2.5, 2.5),
             (-2.5, -1.5), (1.5, 2.5), (-2.5, -.5), (.5, 2.5), (-2.5, 2.5), (-2.5, 2.5), (-2.5, 2.5), (-2.5, -1.5),
             (-0.5, 0.5), (1.5, 2.5), ]
    box_y = [(-2.5, -1.5), (-0.5, 0.5), (1.5, 2.5), (-2.5, 2.5), (-2.5, 1.5), (-2.5, 1.5), (-2.5, 2.5), (-2.5, -1.5),
             (-0.5, 0.5), (1.5, 2.5),
             (2.5, 2.5), (2.5, 2.5), (1.5, 1.5), (1.5, 1.5), (0.5, 0.5), (-0.5, -0.5), (-1.5, -1.5), (-2.5, -2.5),
             (-2.5, -2.5), (-2.5, -2.5), ]
    fields = []

    for (b_x_i, b_y_i) in zip(box_x, box_y):  # Generate mesh refinement boxes around the absortion boundaries
        field = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(field, "VIn", lc / sigma_t)
        gmsh.model.mesh.field.setNumber(field, "VOut", lc)
        gmsh.model.mesh.field.setNumber(field, "XMin", b_x_i[0] - refinement_box_width)
        gmsh.model.mesh.field.setNumber(field, "XMax", b_x_i[1] + refinement_box_width)
        gmsh.model.mesh.field.setNumber(field, "YMin", b_y_i[0] - refinement_box_width)
        gmsh.model.mesh.field.setNumber(field, "YMax", b_y_i[1] + refinement_box_width)
        gmsh.model.mesh.field.setNumber(field, "Thickness", transition_thickness)
        fields.append(field)
    # Finally, let's use the minimum of all the fields as the background mesh field:
    min_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(min_field)
    gmsh.model.occ.synchronize()

    # mark boundaries
    boundary_lines = gmsh.model.occ.getEntities(1)  # 1 specifies edges (lines)
    surface = gmsh.model.occ.getEntities(2)  # 1 specifies edges (lines)

    boundary_marker = gmsh.model.addPhysicalGroup(1, [l[1] for l in boundary_lines[:4]])
    gmsh.model.setPhysicalName(1, boundary_marker, "void")
    surface_marker = gmsh.model.addPhysicalGroup(2, [s[1] for s in surface])
    gmsh.model.setPhysicalName(2, surface_marker, "domain")
    gmsh.model.occ.synchronize()

    gmsh.model.mesh.generate(2)

    # gmsh.write(filename + ".msh")
    # os.system('gmsh ' + filename + '.msh -2 -format su2 -save_all')

    gmsh.write(filename + ".vtk")
    gmsh.write(filename + ".su2")
    # gmsh.write(filename + ".geo_unrolled")
    # Launch the GUI to see the results:
    # if '-nopopup' not in sys.argv:
    #    gmsh.fltk.run()
    # os.system('gmsh ' + filename + '.geo_unrolled -2 -format su2 -save_all')

    gmsh.finalize()  # finalize gmsh

    return 0


if __name__ == '__main__':
    print("---------- Start Creating the mesh ------------")
    print("Parsing options")
    # --- parse options ---
    parser = OptionParser()
    parser.add_option("-o", "--output_name", dest="output_name", default="../lattice_refined")
    parser.add_option("-c", "--char_length", dest="char_length", default=0.5)
    parser.add_option("-s", "--sigma_t", dest="sigma_t", default=11)
    (options, args) = parser.parse_args()
    options.output_name = str(options.output_name)
    options.char_length = float(options.char_length)
    options.sigma_t = float(options.sigma_t)
    main(options.output_name, options.char_length, options.sigma_t)
