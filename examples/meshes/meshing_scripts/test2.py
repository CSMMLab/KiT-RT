import numpy
import meshio
import gmsh
import pygmsh

resolution = 0.01
# Channel parameters
L = 2.2
H = 0.41
c = [0.2, 0.2, 0]
r = 0.05

domain_edge_pts = [[-3.5, -3.5, 0], [3.5, -3.5, 0], [3.5, 3.5, 0], [-3.5, 3.5, 0]]

# Initialize empty geometry using the build in kernel in GMSH
geometry = pygmsh.geo.Geometry()
# Fetch model we would like to add data to
model = geometry.__enter__()
# Add circle
# circle = model.add_circle(c, r, mesh_size=resolution)

# Add points that bound the domain
points = [model.add_point(domain_edge_pts[0], mesh_size=resolution),
          model.add_point(domain_edge_pts[1], mesh_size=resolution),
          model.add_point(domain_edge_pts[2], mesh_size=resolution),
          model.add_point(domain_edge_pts[3], mesh_size=resolution)]
# Add lines between all points creating the rectangle
channel_lines = [model.add_line(points[i], points[i + 1])
                 for i in range(-1, len(points) - 1)]
# Create a line loop and plane surface for meshing
channel_loop = model.add_curve_loop(channel_lines)
plane_surface = model.add_plane_surface(channel_loop)

# Call gmsh kernel before add physical entities
model.synchronize()

volume_marker = 6
model.add_physical([plane_surface], "Volume")
model.add_physical([channel_lines[0], channel_lines[1], channel_lines[2], channel_lines[3]], "void")

geometry.generate_mesh(dim=2)
gmsh.write("mesh_test2.msh")
gmsh.write("../mesh_test2.su2")
gmsh.write("mesh_test2.vtk")

gmsh.clear()
geometry.__exit__()
