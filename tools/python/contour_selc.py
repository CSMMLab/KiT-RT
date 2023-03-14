#Importing libraries
import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv
from PIL import Image
from skimage import io, filters,feature
from skimage.measure import find_contours, approximate_polygon, \
    subdivide_polygon

#Funtion to retrieve the indices of contours between the selected input points.
#Even number of points are selected on the plot using gplot and then using a simple for and if loop
#all the contours within the box created by two successive points. 
def contour_idx(contours, x): 

	""" 
	As we want to use selected contours for generating the mesh, contour_idx extracts the 
	indices of all those contours that lie within a box created by successive points in the 
	set x. 

	Parameters
	----------
	contours : list 
			list of contours, where one contour is an array of size (n,2)
	x : array
			array containing coordinates of points used in creating boxes to search for contours within
	
	Returns
	-------
	out: list
		list containing the indices of selected contours

	"""



	x_num = int(len(x)/2) #half the number of points in the list x
	contour_index = [] #Creating empty list to store indices of the boxed-in contours
	for k  in range(len(contours)):
		c = contours[k][:,0,:] #Extracting all the points in a contour
		nk = len(c) #number of points in a given contour
		for i in range(nk):
		    ckij = c[i] #i^th point in the k^th contour
		    for p in range(0,x_num):
			    if x[2*p][0] < ckij[0] and ckij[0] < x[2*p+1][0]:#Checking if the x-coordinate of the point is within the selected point 
			        if x[2*p][1] > ckij[1] and ckij[1] > x[2*p+1][1]: #Checking if the y-coordinate of the point is within the selected point
			            contour_index.append(k) #appending the contour index k if one of the points of the contour lies within the box
			    else:
			        continue
	return contour_index





	

def contour_selc(img_file, method, num_points, sigma = 2, threshold = 0, max_val = 255):

	""" 
	To be able to contorl the mesh refinement in different regions of an image we need to identify different 
	regions in the image. The function offers two different contouring, watershed seperation or canny contouring
	from OpenCV. 

	Parameters
	----------

	img_file	: jpg or png file
				Image file that is to be used to create a mesh
	method		: str, 'watershed' or 'canny'
	num_points	: int, number of points (even) to be used for selecting contours
	sigma		: float, default = 2
				The level of details wanted in the canny edge detection algorithm
	threshold	: int, default = 0
				The minimum threshold grayscale value used in creating a binary image file
	max_val 	: int, default = 255
				maximum value of grayscale value to be used for creating a binary file 

	
	Returns
	-------
	out:
		contour_index	: list, 
						indices of the contours selected
		contours 		: list of array,
						list of all the contours (as numpy array) in the image according to method
		hierarchy		: numpy array,
						hierarchy tree of the different contours 
		dim 			: list,
						dimensions of the image
 
	"""

	if method == 'watershed':
		img = cv.imread(img_file)
		b,g,r = cv.split(img)
		rgb_img = cv.merge([r,g,b])

		gray = cv.cvtColor(img,cv.COLOR_BGR2GRAY)
		ret, thresh = cv.threshold(gray,0,255,cv.THRESH_BINARY_INV+cv.THRESH_OTSU)

		# noise removal
		kernel = np.ones((2,2),np.uint8)
		#opening = cv.morphologyEx(thresh,cv.MORPH_OPEN,kernel, iterations = 2)
		closing = cv.morphologyEx(thresh,cv.MORPH_CLOSE,kernel, iterations = 2)

		# sure background area
		sure_bg = cv.dilate(closing,kernel,iterations=3)

		# Finding sure foreground area
		dist_transform = cv.distanceTransform(sure_bg,cv.DIST_L2,3)

		# Threshold
		ret, sure_fg = cv.threshold(dist_transform,0.1*dist_transform.max(),255,0)

		# Finding unknown region
		sure_fg = np.uint8(sure_fg)
		unknown = cv.subtract(sure_bg,sure_fg)

		#dimenson of the image
		dim = np.shape(unknown)

		_, binary = cv.threshold(unknown, threshold, max_val, cv.THRESH_BINARY_INV)
		#plt.imshow(binary, cmap="gray")
		binary = np.uint8(binary)

		#image = cv.imread("G:\\HiWi\\SLphantom.jpg")

		contours, hierarchy = cv.findContours(binary, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)

		plt.figure(figsize=(10,10))

		plt.imshow(binary)
		#n = 2*input('Enter the number of points to select the contours to include',)
		x = plt.ginput(num_points)


		#list containing indices of the contours between the selected points
		contour_index = contour_idx(contours, x)

	if method == 'canny':
		img = Image.open(img_file)
		#img = Image.open('G:\\HiWi\\SLphantom.png')
		imgarray = img.convert('LA')
		imgmat = np.array(imgarray.getdata(band=0))
		imgmat.shape = (imgarray.size[1],imgarray.size[0])
		imgmat_norm = imgmat/255

		edges = feature.canny(imgmat_norm,sigma = sigma)
		m,n = np.shape(edges)
		edges_num = np.zeros((m,n))
		for i in range(m):
		    for j in range(n):
		        if edges[i,j] == True:
		            edges_num[i,j] = 1 
		            edges_num[i-1,j] = 1

		#dimenson of the image
		dim = np.shape(edges_num)

		_, binary = cv.threshold(edges_num, threshold, max_val, cv.THRESH_BINARY_INV)
		#plt.imshow(binary, cmap="gray")
		binary = np.uint8(binary)

		image = cv.imread(img_file)
		#image = cv.imread("G:\\HiWi\\SLphantom.jpg")

		contours, hierarchy = cv.findContours(binary, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)

		plt.figure(figsize=(10,10))
		image = cv.drawContours(image, contours, -1, (0,255, 0), 2)

		plt.imshow(binary)
		#n = 2*input('Enter the number of points to select the contours to include',)
		x = plt.ginput(num_points)


		contour_index = contour_idx(contours, x)

	return np.unique(contour_index),contours, hierarchy, dim

# Function to create an list for each contour containing indices of the  
# contours contained within it. For example, if we have a four circles indexed 1,2,3 and 4
# such that 1 contains 2, 3 and 4. circles 2 and 3 exclusive and 2 contains circle 4 within it.
# Then the function returns [[2,3,4],[4],[],[]] where the last two lists are empty as circles 3 and 4 do not
# contain any other circle within it. We can also interpret these relations as a tree.


def contour_tree(hierarchy):

	"""
	The function takes a hierarchy tee as an input and returns, for each contour, the list of contours that are
	contained in it.

	Parameters
	----------

	hierarchy	: array,
				An array containing the (next node, previous node, parent node, child node) information for each
				node

	Return
	------

	out: list
		a list containing the list of contours contained within a particular contour

	"""



	hirchy = hierarchy[0,:,:]
	m,n = np.shape(hirchy)
	holes_list = [[] for k in range(m)]
	for i in range(m):
		h = hirchy[i,:]
		if h[0] == -1 and h[1] == -1:
			holes_list[i].append(h[3])
			holes_list[h[2]].append(i)
		else:
			if h[3] == -1:
				if h[2] == -1:
					if h[1] != -1:
						holes_list[i].append(h[1])
						holes_list[h[1]].append(i)
				else:
					holes_list[i].append(h[2])
			else:
				holes_list[h[3]].append(i)
	return holes_list

