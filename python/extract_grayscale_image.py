# importing PIL
from PIL import Image
import pydicom
import numpy as np
from pydicom.pixel_data_handlers.util import apply_modality_lut
from copy import deepcopy

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
        img = Image.open(image_name).convert('L') #image data
        I = np.asarray(img) # image as greyscale
        I = I/255; # rescale values to [0,1] (1.85 is the density of bone, but it is set in problem::Radioct)
        J = deepcopy(np.flipud(I))
        dimensions = (4,4) # [cm]
        return J , dimensions

