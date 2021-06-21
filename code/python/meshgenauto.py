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
from contour_selc import contour_selc, contour_structure

def Mesh(contours, cont_idx,mesh_holes, dim, indices, lcar, remove_low_dim_cells = True):
	with pg.geo.Geometry() as geom:
		m,n = dim #dimensions of the CT scan image
		
		#creating boundaray point elements in geom
		b1 = geom.add_point([0,0],20) 
		b2 = geom.add_point([0,m],20) 
		b3 = geom.add_point([n,0],20) 
		b4 = geom.add_point([n,m],20)
		
		#using the boundary point elements to create line elements corresponding to image boundary
		line1 = geom.add_line(b1,b2)
		line2 = geom.add_line(b2,b4)
		line3 = geom.add_line(b4,b3)
		line4 = geom.add_line(b3,b1)
		
		#we make a loop element in geom which makes the boundary of the CT scan
		lines = geom.add_curve_loop([line1,line2,line3,line4])
		
		#creating loop elements in geom corresponding to each contour
		loops = [] #empty list to which we append the loop elements created using the contour points
		for k in range(len(cont_idx)):
			l = []	#we append the point elements of the k^th contour to l
			index = indices[k] #we use a subsection of the points of the contour
			cont = contours[cont_idx[k]][:,0,:] #copying the contour points in the variable cont	
			#print(k,len(cont),len(index))
			for i in range(len(index)):
				p = geom.add_point(cont[index[i]],lcar[k]) #creating the point element in geom of the k^th contour
				l.append(p)
				
			#creating a spline element in geom using the points in the list l. For creating a curve_loop element it is required that
			#the first and the last point elements of the two splines must match not only in value but id as well.
			#For example, let say we have points x1,x2,x3,x4,x5,x6,x7 and x8 point elements in geom such that x4 & x5 and x1& x8 are point elements created
			#using the same coordinates (a,b) but have a different id in geom. We create three spline elements 
			#s1 = geom.add_bspline([x1,x2,x3,x4]), s2 = geom.add_bspline([x5,x6,x7,x8]) and s3 = geom.bspline([x4,x6,x7,x1]).
			#Then geom.add_curveloop([s1,s2]) will create an assertion error and instead we must use geom.add_curve_loop([s1,s3])
			s1_l = geom.add_bspline(l) 
			s2_l = geom.add_spline([l[-1],l[0]])

			l_loop = geom.add_curve_loop([s1_l,s2_l])
			loops.append(l_loop)
		
		#creating the boundary surface element in geom with all contours as holes. This surface element can be meshed
		geom.add_plane_surface(lines, holes = loops)	

		#We create surface elements in the mesh for each loop (k). Using the hierarchy structure of the loops we 
		#specify the loops (i1, i2,...,ir) that are contained the loop k and they form holes in surface element 
		#corresponding to loop k

		for k in range(len(cont_idx)):
			holes_idx = mesh_holes[cont_idx[k]] #loading the indices of the loops that form holes in the kth loop
			holes = [] #empty list to append the loop elements according to the indices
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
				
				
		mesh = geom.generate_mesh() #generate the mesh element 
		
		if remove_low_dim_cells == True:
			mesh.remove_orphaned_nodes()
			mesh.remove_lower_dimensional_cells()

		return mesh


def mesh_from_image(image_path, mesh_path,num_points, method = 'watershed', small_contour = False, default_lcar = True, plot_contours = True):

	cont_idx,contours,hierarchy, dim = contour_selc(image_path,method,num_points)

	mesh_holes = contour_tree(hierarchy)

	if small_contour == False:
		del_list = []
		for k in range(len(cont_idx)):
			if len(contours[cont_idx[k]]) < 10:
				del_list.append(k)

		cont_idx = np.delete(cont_idx,del_list)

	if default_lcar == True:
		indices = []
		lcar = []
		for i in range(len(cont_idx)):
			if len(contours[cont_idx[i]]) <= 25:
				index = np.arange(0,len(contours[cont_idx[i]]),1)
				lcar_val = 10
			elif 25 < len(contours[cont_idx[i]]) <= 50:
				index = np.arange(0,len(contours[cont_idx[i]]),2)
				lcar_val = 10
			elif 50 < len(contours[cont_idx[i]]) <= 200:
				index = np.arange(0,len(contours[cont_idx[i]]),5)
				lcar_val = 20
			elif 200 < len(contours[cont_idx[i]]) <= 1000:
				index = np.arange(0,len(contours[cont_idx[i]]),10)
				lcar_val = 25
			else:
				index = np.arange(0,len(contours[cont_idx[i]]),15)
				lcar_val = 40
			indices.append(index)
			lcar.append(lcar_val)

	for i in range(len(cont_idx)):
		plt.scatter(contours[cont_idx[i]][:,0,0],contours[cont_idx[i]][:,0,1], label = '{}'.format(cont_idx[i]))
	plt.legend()
	plt.show()

	mesh = Mesh(contours,cont_idx,mesh_holes,dim, indices,lcar,remove_low_dim_cells = True)
	#help(mesh.write)

	mesh.write(mesh_path)


