rm lattice_n20.su2 lattice_n20.con lattice_n40.su2 lattice_n40.con lattice_n80.su2 lattice_n80.con

gmsh lattice_rectangular_n20.geo -2 -format su2 -save_all -o lattice_n20.su2
gmsh lattice_rectangular_n40.geo -2 -format su2 -save_all -o lattice_n10.su2
gmsh lattice_rectangular_n80.geo -2 -format su2 -save_all -o lattice_n10.su2
