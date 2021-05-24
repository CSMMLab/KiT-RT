#Importing libraries
import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv
from PIL import Image
from skimage import io, filters,feature
from skimage.measure import find_contours, approximate_polygon, \
    subdivide_polygon


def contour_selc(img_file,method, num_points,sigma,threshold,max_val):
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

		_, binary = cv.threshold(edges_num, threshold, max_val, cv.THRESH_BINARY_INV)
		#plt.imshow(binary, cmap="gray")
		binary = np.uint8(binary)

		#image = cv.imread("G:\\HiWi\\SLphantom.jpg")

		contours, hierarchy = cv.findContours(binary, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)

		plt.figure(figsize=(10,10))

		plt.imshow(binary)
		#n = 2*input('Enter the number of points to select the contours to include',)
		x = plt.ginput(num_points)


		contour_index = []
		for k  in range(len(contours)):
			c = contours[k][:,0,:]
			nk = len(c)
			for i in range(nk):
				ckij = c[i]
				for p in range(1,num_points):
					if x[p-1][0] < ckij[0] and ckij[0] < x[p][0]:
						if x[p-1][1] > ckij[1] and ckij[1] > x[p][1]:
							contour_index.append(k)
					else:
						continue
	if method == 'canny':
		img = Image.open(img_file) 
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
					edges_num[i-1,j] = 1 #We need to have some thickness to get contours


		_, binary = cv.threshold(edges_num, threshold, max_val, cv.THRESH_BINARY_INV)
		#plt.imshow(binary, cmap="gray")
		binary = np.uint8(binary) #numbers need to be in 8bit

		image = cv.imread(img_file)

		contours, hierarchy = cv.findContours(binary, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)

		plt.figure(figsize=(10,10))
		image = cv.drawContours(image, contours, -1, (0,255, 0), 2)

		plt.imshow(binary)
		#n = 2*input('Enter the number of points to select the contours to include',)
		x = plt.ginput(num_points)


		contour_index = []
		for k  in range(len(contours)):
			c = contours[k][:,0,:]
			nk = len(c)
			for i in range(nk):
				ckij = c[i]
				for p in range(1,num_points):
					if x[p-1][0] < ckij[0] and ckij[0] < x[p][0]:
						if x[p-1][1] > ckij[1] and ckij[1] > x[p][1]:
							contour_index.append(k)
					else:
						continue
		return np.unique(contour_index),contours, hierarchy, dim

# We create  function called 'mesh_holes' to collect the hierarchies of contours
# using the hierarchy tree. We create parent-chid relationships between the different
# contours.

def mesh_holes(hierarchy):
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


#'G:\\HiWi\\Liver_CT.png'
contour_index,contours = contour_selc("G:\\HiWi\\Lung_CT.png",10,2,0,255)
print(contour_index)
for i in range(len(contour_index)):
	plt.scatter(contours[contour_index[i]][:,0,0],contours[contour_index[i]][:,0,1], label = '{}'.format(contour_index[i]))
plt.legend()
plt.show()
