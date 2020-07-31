'''
Author: Steffen Schotthöfer
Description: This is the first prototype of a network to use in the entropy optimizer
https://www.tensorflow.org/tutorials/keras/save_and_load
'''

import numpy as np                   # advanced math library
import random                        # for generating random numbers
import tensorflow as tf


# Define the NN model
# Build the network:
def create_model():
  model = tf.keras.models.Sequential([
    keras.layers.Dense(512, activation='relu', input_shape=(4,)),
    keras.layers.Dropout(0.2),
    keras.layers.Dense(512, activation='relu', input_shape=(512,)),
    keras.layers.Dropout(0.2),
    keras.layers.Dense(4)
  ])

  model.summary()
  model.compile(optimizer='adam',
                loss=tf.losses.MeanSquaredError(),
                metrics=['accuracy'])

  return model


#Prepare the Data
def prepData(Data):
    divide
    into
    input and output
    data


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

    return (xDataTrain,yDataTrain,xDataTest,yDataTest)
