'''
Author: Steffen Schotthöfer
Description: This file is used to explore load a saved pre-trained keras network in tf2.0 format and employ it
'''

# Read test data
# Load the Pandas libraries with alias 'pd'
import pandas as pd
import random
import numpy as np
import os
import tensorflow as tf
from tensorflow import keras

# get some data
dataFrameInput = pd.read_csv("trainNN.csv")
# Preview the first 5 lines of the loaded data
dataFrameInput.head()

# divide into input and output data

idx = 0
xDataList = list()
yDataList = list()
data = dataFrameInput.values  # numpy array

## Somehow t = 0 has an odd number of elements, just cut it out

for row in data:
    if (row[2] > 0):
        if (idx % 2 == 0):
            xDataList.append(row)
        else:
            yDataList.append(row)

        idx = idx + 1

# merge the lists
DataList = list()
for rowX, rowY in zip(xDataList, yDataList):
    DataList.append([rowX, rowY])

# Shuffle data

random.shuffle(DataList)

DataArray = np.asarray(DataList)
print(DataArray.shape)

# Strip off header information, i.e. the first 3 cols
DataArraySlim = DataArray[:, :, 3:]
print(DataArraySlim.shape)

# split in train and test data (ratio 2:1)
DataTrain = DataArraySlim[:2 * int(DataArraySlim.shape[0] / 3)]
DataTest = DataArraySlim[2 * int(DataArraySlim.shape[0] / 3):]

# Split in x (input) and y (output) data
xDataTrain = DataTrain[:, 0, :]
yDataTrain = DataTrain[:, 1, :]
xDataTest = DataTest[:, 0, :]
yDataTest = DataTest[:, 1, :]

print(yDataTrain.shape)
print(yDataTest.shape)



# Load model
new_model = tf.keras.models.load_model('saved_model/my_model')

# Check its architecture
new_model.summary()

# Evaluate the restored model
loss, acc = new_model.evaluate(xDataTest,  yDataTest, verbose=2)
print('Restored model, accuracy: {:5.2f}%'.format(100*acc))

print(new_model.predict(xDataTest).shape)
