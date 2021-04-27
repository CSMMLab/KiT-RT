#Importing libraries
import numpy as np
import matplotlib.pyplot as plt
import cv2 as cv
from PIL import Image
from skimage import io, filters,feature
from skimage.measure import find_contours, approximate_polygon, \
    subdivide_polygon


def contour_selc(img_file, num_points,sigma,thresh,max_val):
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


	_, binary = cv.threshold(edges_num, thresh, max_val, cv.THRESH_BINARY_INV)
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


	print(x)

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
	return np.unique(contour_index),contours


#'G:\\HiWi\\Liver_CT.png'
contour_index,contours = contour_selc("G:\\HiWi\\Lung_CT.png",10,2,0,255)
print(contour_index)
for i in range(len(contour_index)):
	plt.scatter(contours[contour_index[i]][:,0,0],contours[contour_index[i]][:,0,1], label = '{}'.format(contour_index[i]))
plt.legend()
plt.show()