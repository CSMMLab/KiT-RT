#rm lattice_n20.su2 lattice_n20.con 
#rm lattice_n40.su2 lattice_n40.con 
#rm lattice_n80.su2 lattice_n80.con
rm quarter_sym_hohlraum_p02.su2 quarter_sym_hohlraum_p02.con
rm quarter_sym_hohlraum_p01.su2 quarter_sym_hohlraum_p01.con
rm quarter_sym_hohlraum_p005.su2 quarter_sym_hohlraum_p005.con
rm quarter_sym_hohlraum_p0025.su2 quarter_sym_hohlraum_p0025.con
rm quarter_sym_hohlraum_p00125.su2 quarter_sym_hohlraum_p00125.con
rm quarter_sym_hohlraum_p00075.su2 quarter_sym_hohlraum_p00075.con


gmsh quarter_sym_hohlraum_p02.geo -2 -format su2 -save_all -o quarter_sym_hohlraum_p02.su2
gmsh quarter_sym_hohlraum_p01.geo -2 -format su2 -save_all -o quarter_sym_hohlraum_p01.su2
gmsh quarter_sym_hohlraum_p005.geo -2 -format su2 -save_all -o quarter_sym_hohlraum_p005.su2
gmsh quarter_sym_hohlraum_p0025.geo -2 -format su2 -save_all -o quarter_sym_hohlraum_p0025.su2
gmsh quarter_sym_hohlraum_p00125.geo -2 -format su2 -save_all -o quarter_sym_hohlraum_p00125.su2
gmsh quarter_sym_hohlraum_p00075.geo -2 -format su2 -save_all -o quarter_sym_hohlraum_p00075.su2
