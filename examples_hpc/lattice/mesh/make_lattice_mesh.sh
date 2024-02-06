#rm lattice_n20.su2 lattice_n20.con 
#rm lattice_n40.su2 lattice_n40.con 
#rm lattice_n80.su2 lattice_n80.con
rm lattice_n160.su2 lattice_n160.con

#gmsh lattice_rectangular_n20.geo -2 -format su2 -save_all -o lattice_n20.su2
#gmsh lattice_rectangular_n40.geo -2 -format su2 -save_all -o lattice_n40.su2
#gmsh lattice_rectangular_n80.geo -2 -format su2 -save_all -o lattice_n80.su2
gmsh lattice_rectangular_n160.geo -2 -format su2 -save_all -o lattice_n160.su2
