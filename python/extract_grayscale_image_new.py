# importing PIL
from PIL import Image
import pydicom
import numpy as np
from pydicom.pixel_data_handlers.util import apply_modality_lut
from copy import deepcopy
import cv2

def extract(image_name):

    ending = image_name.split('.')[-1];

    if(ending == 'dcm'): #dicom ct image
        dataset = pydicom.dcmread(image_name)
        hu = apply_modality_lut(dataset.pixel_array, dataset)

        density = (dataset.RescaleIntercept + dataset.RescaleSlope*hu)/1000 +1
        (spacingX, spacingY) = dataset.PixelSpacing
        (sizeX,sizeY) = density.shape
        dimensions = (spacingX*sizeX, spacingY*sizeY)
        return density, dimensions

    else: # jpg or png
        img = cv2.imread(image_name, cv2.IMREAD_UNCHANGED)
        scale_percent = 15
        print(img.shape)
        width = int(img.shape[1] * scale_percent / 100)
        height = int(img.shape[0] * scale_percent / 100)
        dim = (width, height)
        img = cv2.resize(img, dim, interpolation = cv2.INTER_AREA)
        print(img.shape)
        I = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY);#np.asarray(img) # image as greyscale
        I = I/255; # rescale values to [0,1] (1.85 is the density of bone, but it is set in problem::Radioct)
        J = deepcopy(np.flipud(I))
        dimensions = (6,6) # [cm]
        return J , dimensions

