import numpy as np
import pygmsh
import os
import pyvista

lc = 100.0  # Characteristic length of mesh
xmin, xmax = 0.0, 2000.0  # X axis boundaries
zmin, zmax = -500.0, -2500.0  # Z axis boundaries

inj_z = -1500.0  # Depth of injection
flt_offset = 50.0  # Offset of fault
flt_thick = 10.0  # Thickness of fault
tana = np.tan(np.deg2rad(80.0))  # Tangeant of dipping angle of fault
dist = 500.0 - 0.5 * flt_thick  # Distance from injection point (0.0, -1500.0) to left wall of fault

bnd_thick = 10.0  # Thickness of boundary elements

bound_right_pts = [
    [xmax, zmin, 0.0],
    [xmax + bnd_thick, zmin, 0.0],
    [xmax + bnd_thick, zmax, 0.0],
    [xmax, zmax, 0.0],
    [xmax, -1700.0 + flt_offset, 0.0],
    [xmax, -1550.0 + flt_offset, 0.0],
    [xmax, -1450.0 + flt_offset, 0.0],
    [xmax, -1300.0 + flt_offset, 0.0],
]
bound_top_pts = [
    [xmin, zmin, 0.0],
    [dist + (zmin - inj_z) / tana, zmin, 0.0],
    [dist + (zmin - inj_z) / tana + flt_thick, zmin, 0.0],
    [xmax, zmin, 0.0],
    [xmax + bnd_thick, zmin, 0.0],
    [xmax + bnd_thick, zmin + bnd_thick, 0.0],
    [dist + (zmin - inj_z + bnd_thick) / tana + flt_thick, zmin + bnd_thick, 0.0],
    [dist + (zmin - inj_z + bnd_thick) / tana, zmin + bnd_thick, 0.0],
    [xmin, zmin + bnd_thick, 0.0],
]
bound_bot_pts = [
    [xmin, zmax, 0.0],
    [dist + (zmax - inj_z) / tana, zmax, 0.0],
    [dist + (zmax - inj_z) / tana + flt_thick, zmax, 0.0],
    [xmax, zmax, 0.0],
    [xmax + bnd_thick, zmax, 0.0],
    [xmax + bnd_thick, zmax - bnd_thick, 0.0],
    [dist + (zmax - inj_z - bnd_thick) / tana + flt_thick, zmax - bnd_thick, 0.0],
    [dist + (zmax - inj_z - bnd_thick) / tana, zmax - bnd_thick, 0.0],
    [xmin, zmax - bnd_thick, 0.0],
]


def add_block(x0, y0, length):
    coords = np.array([
        [x0, y0, 0.0],
        [x0 + length, y0, 0.0],
        [x0 + length, y0 + length, 0.0],
        [x0, y0 + length, 0.0]
    ])
    return coords


lc = 0.1
x_inner = [-2.5, -2.5, -2.5, -1.5, -1.5, -0.5, -0.5, .5, .5, 1.5, 1.5, 1.5]
y_inner = [1.5, -0.5, -2.5, 0.5, -1.5, -0.5, -2.5, 0.5, -1.5, 1.5, -0.5, -2.5]

boxes = []
for x_i, y_i in zip(x_inner, y_inner):
    box = add_block(x_i, y_i, 1)
    boxes.append(box)

domain_edge_pts = [[-3.5, -3.5, 0], [3.5, -3.5, 0], [3.5, 3.5, 0], [-3.5, 3.5, 0]]
with pygmsh.geo.Geometry() as geo:
    # Define polygons

    domain = geo.add_polygon(domain_edge_pts, mesh_size=lc)
    # geo.boolean_fragments([domain], [])

    # geo.add_physical(domain, "domain")
    geo.add_physical(domain.lines, label="void")

    # Remove duplicate entities
    mesh = geo.generate_mesh(dim=2)  # algorithm=6

    # Export the mesh for post-processing
    mesh.write("mesh.vtk")
    mesh.write("../mesh2.su2")
    mesh.write("../mesh2.geo_unrolled")
    # mesh_code = geo.get_code()
    # with open("mesh.geo", "w") as mesh_file:
    #    mesh_file.write(mesh_code)
    # os.system('gmsh ' + 'mesh.geo -2 -format su2 -save_all')

pyvista.set_plot_theme("document")

p = pyvista.Plotter(window_size=(800, 800))
p.add_mesh(
    mesh=pyvista.from_meshio(mesh),
    scalar_bar_args={"title": "Materials"},
    show_scalar_bar=True,
    show_edges=True,
)
p.view_xy()
p.show()
