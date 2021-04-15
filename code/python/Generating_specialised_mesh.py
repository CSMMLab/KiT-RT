#!/usr/bin/env python
# coding: utf-8

# # Specific Mesh for SLphantom image

# In[ ]:


#Importing necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pygmsh as pg
from PIL import Image
from scipy import ndimage as ndi

from skimage import io, filters,feature
from skimage.measure import find_contours, approximate_polygon,     subdivide_polygon
import meshio
import cv2 as cv
import pymesh
import matplotlib.colors as mcolors


# In[ ]:


#Loading the image file
img = Image.open('G:\\HiWi\\SLphantom.png')
imgarray = img.convert('LA')
imgmat = np.array(imgarray.getdata(band=0))
imgmat.shape = (imgarray.size[1],imgarray.size[0])
#normalising the image pixel values to be between 0 and 1
imgmat_norm = imgmat/255 


# We use the OpenCV's python extension librarz cv2 library for edge detection and finding contours of the image (https://docs.opencv.org/4.5.1/d6/d00/tutorial_py_root.html). The idea is to convert a grey scale image into a binary matrix using cv.threshold(image, thresh, maxval, \*args) and setting a particular thresh value. We then feed this binary file into cv.findContours(binary,\*args) function to get contours in the binary file as well as the hierarchy of the countours. 

# In[ ]:


image = cv.imread("G:\\HiWi\\SLphantom.png")

#Creating grey scale image
image = cv.cvtColor(image, cv.COLOR_BGR2RGB)
gray = cv.cvtColor(image, cv.COLOR_RGB2GRAY)

#Creating binary image file using threshold 
thresh = 50
_, binary = cv.threshold(gray, thresh, 255, cv.THRESH_BINARY_INV)
#plt.imshow(binary, cmap="gray")
#plt.show()

#Finding countours in the binary image as well as the hierarchy of the countours
contours, hierarchy = cv.findContours(binary, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)

#Plotting the countours on the original image
#image = cv.drawContours(image, contours, -1, (0, 255, 0), 2)
#plt.figure(figsize=(10,10))
#plt.imshow(image)
#plt.show()


# ### Limiations
# 1. As the contours that we find are dependant on the threshold value that we use, we get different objects in the image on the for different threshold values 
# For example, for thresh = 50 we get,
# ![title](img/binary_image_SLp.png)
# So we set different threshold values (50,75,85,250) to retireve all the different contours in the image.
# 
# 2. This method still does not give certain contours in the image and leaves out certain portions of objects that intersect.

# In[ ]:


#from threshhold value 250
layer1 = contours[1][:,0,:]
layer2 = contours[2][:,0,:]


# In[ ]:


#threshhold value = 50
layer3 = contours[2][:,0,:]
layer4 = contours[3][:,0,:]


# In[ ]:


#threshhold = 75
layer5 = contours[3][:,0,:]
layer6 = contours[4][:,0,:]
layer7 = contours[5][:,0,:]
layer8 = contours[6][:,0,:]
layer9 = contours[7][:,0,:]


# In[ ]:


#threshhold=85
layer10 = contours[3][:,0,:]


# After using different threshold values and extracting different contours, we get the following plot
# ![text](img/extract_contours.png)

# Now we use the python interface of the gmsh (pygmsh) to create the mesh. We need to feed in the in the form of curved loops and then mesh them. For this we follow the following procedure:
# 1. create a pg.geo.Geometry() instance as geom
# 2. For including a particular we add each point from the list by iterating over an index list with a parameter lcar which would determine the mesh size
#     1. We do not want to use all the points on the contour as that will restrict our freedom in choosing the size of the mesh. So we iterate over a index list which is a subset of the entire list
# 3. Use these points to create a spline or bspline
# 4. Using the splines create a curved loop instance in geom. This curved loop can now be meshed by adding it as a surface. Note: While creating the loop we must ensure that the end points of the two splines must be the exact same, which means that they must have the same object id as well.
# 5. Each surface that we include into geom using geom.add_plane_surface(curvedloop,holes = [ ])
# 6. We have to make further adjustments when two surfaces intersect for instance l9 & l4 and l8 & l4
#     1. identify the common edge
#     2. create a spline of the common edge similar to any other contour
#     3. keeping in mind the end point assertion create a revised loop and replace the older loops
# 7. create plane surfaces that can be meshed by adding appropriate holes. For example, the surface l1 has curved loop l2 as a hole and surface l2 has l3-l10 loops as holes. 
# 

# In[ ]:


#indexing sets for each of the contours
index1 = np.arange(0,len(layer1),28)
index2 = np.arange(0,len(layer2),28)
index3 = np.arange(0,len(layer3),10)
index4 = np.arange(0,len(layer4_re1),10)
index5 = np.arange(0,len(layer5),2)
index6 = np.arange(0,len(layer6),3)
index7 = np.arange(0,len(layer7),3)
index8 = np.arange(0,len(layer8_re),3)
index9 = np.arange(0,len(layer9_re),10)
index10 = np.arange(0,len(layer10),2)
index_m1 = np.arange(0,len(mid1),4)
index_m2 = np.arange(0,len(mid2),2)
index_l4mid = np.arange(0,len(layer4_mid),4)


# In[ ]:


with pg.geo.Geometry() as geom:
    lcar = 20
    lcar2 = 25
    lcar1 = 75
    lcar3 = 20
    lcar4 = 10
    lcar5 = 15
    lcar6 = 15
    lcar7 = 15
    lcar8 = 5
    lcar9 = 50
    lcar10 = 10
    
    #creating the outer boundary of the mesh
    b1 = geom.add_point([0,0],lcar) #adding the corner points as pygmsh point instance
    b2 = geom.add_point([0,1200],lcar)
    b3 = geom.add_point([1200,0],lcar)
    b4 = geom.add_point([1200,1200],lcar)
    
    line1 = geom.add_line(b1,b2)  #Creating lines using point instances of the four corners
    line2 = geom.add_line(b2,b4)
    line3 = geom.add_line(b4,b3)
    line4 = geom.add_line(b3,b1)
    lines = geom.add_curve_loop([line1,line2,line3,line4]) #creating a curved loop using the 4 lines
    
    #Creating the outermost contour of the image
    l1 = [] # creating list to store the points instances created in geom of the contour
    for i in range(len(index1)):
        p = geom.add_point(layer1[index1[i]],lcar) #creating the point instance in geom
        l1.append(p) 
    
    s1_l1 = geom.add_bspline(l1) #adding a bspline instance in geom that gows through all the points
    s2_l1 = geom.add_spline([l1[-1],l1[0]]) #creating a second spline instance to use to close the curved loop
    
    l1_loop = geom.add_curve_loop([s1_l1,s2_l1]) #creating a curved loop
    
    #loop 2
    l2 = []
    for i in range(len(index2)):
        p = geom.add_point(layer2[index2[i]],lcar2)
        l2.append(p)
    
    s1_l2 = geom.add_bspline(l2)
    s2_l2 = geom.add_spline([l2[-1],l2[0]])
    
    l2_loop = geom.add_curve_loop([s1_l2,s2_l2])
    
    geom.add_plane_surface(l1_loop,holes = [l2_loop])
   
    #loop 3
    l3 = []
    for i in range(len(index3)):
        p = geom.add_point(layer3[index3[i]],lcar3)
        l3.append(p)
    
    s1_l3 = geom.add_bspline(l3)
    s2_l3 = geom.add_spline([l3[-1],l3[0]])
    
    l3_loop = geom.add_curve_loop([s1_l3,s2_l3])
    
    #intersection between l4 and l9
    m1 = []
    for i in range(len(index_m1)):
        p = geom.add_point(mid1[index_m1[i]],lcar4)
        m1.append(p)
    
    s_m1 = geom.add_bspline(m1)
    
    #intersection between l4 and l8
    m2 = []
    for i in range(len(index_m2)):
        p = geom.add_point(mid2[index_m2[i]],lcar4)
        m2.append(p)
    
    s_m2 = geom.add_bspline(m2)
    
    #the part of l4 between s_m1 and s_m2
    l4_mid = []
    for i in range(len(index_l4mid)):
        p = geom.add_point(layer4_mid[index_l4mid[i]],lcar4)
    
    l4_mid.insert(0,m2[-1])
    l4_mid.insert(len(l4_mid),m1[0])
    s_l4mid = geom.add_bspline(l4_mid)
    
    #loop 4
    l4 = []
    for i in range(len(index4)):
        p = geom.add_point(layer4_re1[index4[i]],lcar4)
        l4.append(p)
    
    l4.insert(0,m1[-1])
    l4.insert(len(l4),m2[0])
    s1_l4 = geom.add_bspline(l4)
    
    l4_loop = geom.add_curve_loop([s1_l4,s_m2,s_l4mid,s_m1])
    
    
    l5 = []
    for i in range(len(index5)):
        p = geom.add_point(layer5[index5[i]],lcar5)
        l5.append(p)
    
    s1_l5 = geom.add_bspline(l5)
    s2_l5 = geom.add_spline([l5[-1],l5[0]])
    
    l5_loop = geom.add_curve_loop([s1_l5,s2_l5])
    
    
    l6 = []
    for i in range(len(index6)):
        p = geom.add_point(layer6[index6[i]],lcar6)
        l6.append(p)
    
    s1_l6 = geom.add_bspline(l6)
    s2_l6 = geom.add_spline([l6[-1],l6[0]])
    
    l6_loop = geom.add_curve_loop([s1_l6,s2_l6])
 
    l7 = []
    for i in range(len(index7)):
        p = geom.add_point(layer7[index7[i]],lcar7)
        l7.append(p)
    
    s1_l7 = geom.add_bspline(l7)
    s2_l7 = geom.add_spline([l7[-1],l7[0]])
    
    l7_loop = geom.add_curve_loop([s1_l7,s2_l7])
    
    l8 = []
    for i in range(len(index8)):
        p = geom.add_point(layer8_re[index8[i]],lcar8)
        l8.append(p)
    
    l8.insert(0,m2[-1])
    l8.insert(len(l8),m2[0])
    s1_l8 = geom.add_bspline(l8)
    
    l8_loop = geom.add_curve_loop([s1_l8,s_m2])
    
    l9 = []
    for i in range(len(index9)):
        p = geom.add_point(layer9_re[index9[i]],lcar9)
        l9.append(p)
    
    l9.insert(0,m1[-1])
    l9.insert(len(l9),m1[0])
    s1_l9 = geom.add_bspline(l9)
    
    l9_loop = geom.add_curve_loop([s1_l9,s_m1])
    
    
    l10 = []
    for i in range(len(index10)):
        p = geom.add_point(layer10[index10[i]],lcar10)
        l10.append(p)
    
    s1_l10 = geom.add_bspline(l10)
    s2_l10 = geom.add_spline([l10[-1],l10[0]])
    
    l10_loop = geom.add_curve_loop([s1_l10,s2_l10])
    
    #creating plane surface instance in geom that can be meshed
    geom.add_plane_surface(lines,holes = [l1_loop]) #accrding to hierarchy we create plane surfaces with holes in them
    
    geom.add_plane_surface(l2_loop,holes = [l3_loop,l4_loop,l5_loop,l6_loop,l7_loop,l8_loop,l9_loop,l10_loop])
    
    geom.add_plane_surface(l3_loop)
    
    geom.add_plane_surface(l4_loop,holes = [l8_loop,l9_loop])
    
    geom.add_plane_surface(l5_loop)
    
    geom.add_plane_surface(l6_loop)
    
    geom.add_plane_surface(l7_loop)
    
    geom.add_plane_surface(l8_loop)
    
    geom.add_plane_surface(l9_loop,holes = [l10_loop])
    
    geom.add_plane_surface(l10_loop)
    
    mesh = geom.generate_mesh() #generating the mesh 


# In[ ]:


mesh.remove_orphaned_nodes() #removing orpaned nodes from the mesh
mesh.remove_lower_dimensional_cells() #removing the curved loops and points (which are lower dimension than triangles) from the mesh
mesh.write("G:\\HiWi\\test1.vtk") #saving the mesh as a .vtk file that can be viewed using ParaView


# ![text](img/mesh_SLphantom_corrected.png)

# In[ ]:




