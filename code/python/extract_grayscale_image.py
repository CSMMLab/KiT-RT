# importing PIL
from PIL import Image
import numpy as np

def extract(image_name):
    img = Image.open(image_name).convert('L') #image data
    I = np.asarray(img) # image as greyscale
    I = I/255; # rescale values to [0,1]
    dimensions = (1,1) # [cm]
    return I , dimensions