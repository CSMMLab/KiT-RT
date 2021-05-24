#Auto mesh generator 
#Importing libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pygmsh as pg
from PIL import Image
from scipy import ndimage as ndi

from skimage import io, filters,feature
from skimage.measure import find_contours, approximate_polygon, \
    subdivide_polygon
import meshio
import cv2 as cv
import pymesh
import matplotlib.colors as mcolors
from contour_selc import contour_selc, mesh_holes

def Mesh(contours, cont_idx,mesh_holes, dim, indices, lcar):
	with pg.geo.Geometry() as geom:
		m,n = dim
		b1 = geom.add_point([0,0],20)
		b2 = geom.add_point([0,m],20)
		b3 = geom.add_point([n,0],20)
		b4 = geom.add_point([n,m],20)
		line1 = geom.add_line(b1,b2)
		line2 = geom.add_line(b2,b4)
		line3 = geom.add_line(b4,b3)
		line4 = geom.add_line(b3,b1)

		lines = geom.add_curve_loop([line1,line2,line3,line4])

		loops = []
		for k in range(len(cont_idx)):
			l = []	
			index = indices[k]
			cont = contours[cont_idx[k]][:,0,:]	
			print(k,len(cont),len(index))
			for i in range(len(index)):
				p = geom.add_point(cont[index[i]],lcar[k])
				l.append(p)

			s1_l = geom.add_bspline(l)
			s2_l = geom.add_spline([l[-1],l[0]])

			l_loop = geom.add_curve_loop([s1_l,s2_l])
			loops.append(l_loop)

		geom.add_plane_surface(lines, holes = loops)	
		#print(loops)

		for k in range(len(cont_idx)):
			holes_idx = mesh_holes[cont_idx[k]]
			holes = []
			for i in range(len(holes_idx)):
				for j in range(len(cont_idx)):
					if holes_idx[i] == cont_idx[j]:
						holes.append(loops[j])
			if not holes:
				geom.add_plane_surface(loops[k])
				print('surface with no holes',cont_idx[k])
			else:
				geom.add_plane_surface(loops[k], holes = holes) 
				print('surface with holes',cont_idx[k])
		mesh = geom.generate_mesh()
		return mesh


cont_idx,contours,hierarchy, dim = contour_selc('G:\\HiWi\\head_CT.jpg','watershed',4,2,0,255)
mesh_holes = mesh_holes(hierarchy)
#mesh_holes[0].remove(59)
del_list = []
for k in range(len(cont_idx)):
	if len(contours[cont_idx[k]]) < 5:
		del_list.append(k)

cont_idx = np.delete(cont_idx,del_list)
		
for i in range(len(cont_idx)):
	plt.scatter(contours[cont_idx[i]][:,0,0],contours[cont_idx[i]][:,0,1], label = '{}'.format(cont_idx[i]))
plt.legend()
plt.show()

indices = []
for i in range(len(cont_idx)):
	if len(contours[cont_idx[i]]) <= 25:
		index = np.arange(0,len(contours[cont_idx[i]]),1)
	elif 25 < len(contours[cont_idx[i]]) <= 50:
		index = np.arange(0,len(contours[cont_idx[i]]),2)
	elif 50 < len(contours[cont_idx[i]]) <= 200:
		index = np.arange(0,len(contours[cont_idx[i]]),5)
	elif 200 < len(contours[cont_idx[i]]) <= 1000:
		index = np.arange(0,len(contours[cont_idx[i]]),10)
	else:
		index = np.arange(0,len(contours[cont_idx[i]]),15)
	indices.append(index)

#lcar = [10,10,5,5,5,5]
lcar = []
for i in range(len(cont_idx)):
	lcar.append(10)

mesh = Mesh(contours,cont_idx,mesh_holes,dim, indices,lcar)
mesh.remove_orphaned_nodes()
#mesh.remove_lower_dimensional_cells()
mesh.write("G:\\HiWi\\autogentest.vtk")

