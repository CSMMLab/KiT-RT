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
#import os
#os.environ.pop("QT_QPA_PLATFORM_PLUGIN_PATH")
import pymesh
import matplotlib.colors as mcolors
from contour_selc import contour_selc, contour_tree
def Mesh(contours, cont_idx,mesh_holes, dim, indices, lcar, remove_low_dim_cells = True):

	"""
	Creates a mesh based on the contours given as input along with containment hierarchy. 

	Parameters
	----------

	contours 	: list,
				list containing a list of contour points
	cont_idx 	: list,
				list of indices of contours to be meshed
	mesh_holes	: list,
				list containing the information regarding the child contours of each contour
	dim 		:list (int,int),
				dimension of the image file
	indices		: list of arrays,
				a list containing the subset of indices of the contours to be used in the mesh
	lcar 		: list,
				the meshing distance parameter
	remove_low_dim_cells : Bool, default True
				Removes the lower dimensional objects from the mesh after meshing

	Returns
	-------
	out: 
		mesh: meshio instance


	"""


	with pg.geo.Geometry() as geom:
		
		m,n = dim #dimensions of the CT scan image

		b1 = geom.add_point([0,0],40) #first boundary point element in the mesh
		b2 = geom.add_point([0,m],40) #second boundary point element in the mesh
		b3 = geom.add_point([n,0],40) #third boundary point element in the mesh
		b4 = geom.add_point([n,m],40) #fourth boundary point element in the mesh
		line1 = geom.add_line(b1,b2) #we create a line element in the mesh between points b1 and b2
		line2 = geom.add_line(b2,b4) #we create a line element in the mesh between points b2 and b4
		line3 = geom.add_line(b4,b3) #we create a line element in the mesh between points b4 and b3
		line4 = geom.add_line(b3,b1) #we create a line element in the mesh between points b3 and b1

		lines = geom.add_curve_loop([line1,line2,line3,line4]) #Creating a loop in the mesh using the four lines created above

		loops = [] #empty list to which we append the loop elements created using the contours
		for k in range(len(cont_idx)):
			l = []	#empty list to append the point elements created using the points in a particular contour
			index = indices[k] #we do not use all the points of a contour but only a portion of it which we obtain using specific indices
			cont = contours[cont_idx[k]][:,0,:]	 #extracting the list of points of a contour
			#print(k,len(cont),len(index)) 
			for i in range(len(index)):
				p = geom.add_point(cont[index[i]],lcar[k]) #creating a point element in the mesh of the k^th contour
				l.append(p) #appending the point element to l 

			s1_l = geom.add_bspline(l) #creating a spline using the points in the list l
			#creating a two point spline element with the first and last point of the contour as the add_curve_loop function requires 
			#that the first and the last points of two spline elements that are to be stitched together must have the same point element
			s2_l = geom.add_spline([l[-1],l[0]]) 

			l_loop = geom.add_curve_loop([s1_l,s2_l]) #creating the loop element of the contour
			loops.append(l_loop) #appending this element to loops

		geom.add_plane_surface(lines, holes = loops) #creating boundary surface element that can now be meshed
		#print(loops)


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
	
def contrast_equ(img_filename,out_filename):
	image = cv.imread(img_filename,0)

	clahe = cv.createCLAHE(clipLimit = 5.0, tileGridSize = (8,8))
	cl1 = clahe.apply(image)
	equ = cv.equalizeHist(cl1)

	cv.imwrite(out_filename,equ)


def mesh_from_image(image_path, mesh_path, num_points, method = 'watershed', small_contour = False, default_lcar = True, plot_contours = True):
	
	"""
	Creates a meshio instance containing the desired mesh based on the given image.

	Parameters
	----------

	image_path	: str,
				location of the image on the device
	mesh_path	: str,
				location where the mesh should be stored with the name of the file. Must be in .vtk format, for example, 
				'C:/Documents/meshautogen.vtk'
	num_points	: int,
				number of points to be used for selecting contours, must be an even number
	method		: str, 
				method to use for contouring the image, default is 'watershed'
	small_contour : Bool,
				If True, all the contours which have lower than 10 points are removed
	default_lcar : Bool or list,
				If True uses the in-built meshing parameters for the mesh, feed this as a list otherwise
	plot_contours : Bool,
				If True plots the selected contours on top of the contoured image, default is False


	Returns
	-------




	"""	

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
	else:
		lcar = default_lcar

	for i in range(len(cont_idx)):
		plt.scatter(contours[cont_idx[i]][:,0,0],contours[cont_idx[i]][:,0,1], label = '{}'.format(cont_idx[i]))
	plt.legend()
	plt.show()

	plt.close()


	mesh = Mesh(contours,cont_idx,mesh_holes,dim, indices,lcar,remove_low_dim_cells = True)
	#help(mesh.write)

	mesh.write(mesh_path)

	
contrast_equ('G:\\HiWi\\Prostate_1.png','G:\\HiWi\\Prostate_1_conteq.png')

mesh_from_image('G:\\HiWi\\Prostate_1_conteq.png','G:\\HiWi\\autogentest.vtk',4)



