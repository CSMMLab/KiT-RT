# importing PIL
from PIL import Image
import pydicom
import numpy as np
from pydicom.pixel_data_handlers.util import apply_modality_lut
from copy import deepcopy

def main():
    image_name = "mini_phantom.png"

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
        I = I/255; # rescale values to [0,1]
        J = deepcopy(np.flipud(I))
        np.savetxt("foo.csv", J, delimiter=",")
        dimensions = (1,1) # [cm]
        return J , dimensions


if __name__ == '__main__':
    main()
