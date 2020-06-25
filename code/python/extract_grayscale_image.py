# importing PIL
from PIL import Image
import numpy as np
from copy import deepcopy

def extract(image_name):
    img = Image.open(image_name).convert('L') #image data
    I = np.asarray(img) # image as greyscale
    I = I/255; # rescale values to [0,1]
    J = deepcopy(np.flipud(I))
    dimensions = (1,1) # [cm]
    return J, dimensions