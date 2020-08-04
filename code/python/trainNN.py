'''
Author: Steffen Schotthöfer
Description: This file contains the Training routine for the entropy NN
'''

# Imports
import pandas as pd
import random                        # for generating random numbers
import numpy as np
import tensorflow as tf
from tensorflow import keras
import matplotlib.pyplot as plt      # MATLAB like plotting routines

def main():
    (xDataTrain,yDataTrain,xDataTest,yDataTest) = prepare_data("trainNN.csv")
    print(yDataTrain.shape)
    print(yDataTest.shape)

    model = create_model()

    cb_list = create_callbacks()
    # Do the training
    history = model.fit(xDataTrain, yDataTrain, validation_split=0.33, epochs=500, batch_size=1000, verbose=1)

    # plot training results
    fig, axs = plt.subplots(2)
    print(history.history.keys())
    axs[0].plot(history.history['accuracy'])
    axs[0].plot(history.history['val_accuracy'])
    axs[0].set_title('model accuracy')
    axs[0].set_ylabel('accuracy')
    axs[0].set_xlabel('epoch')
    axs[0].legend(['train_acc', 'val_acc'], loc='upper left')
    #axs[0].show()
    # summarize history for loss
    axs[1].plot(history.history['loss'])
    axs[1].plot(history.history['val_loss'])
    axs[1].set_title('model loss')
    axs[1].set_ylabel('loss')
    axs[1].set_xlabel('epoch')
    axs[1].legend(['train_loss', 'val_loss'], loc='upper left')
    #axs[1].show()
    plt.show()
    # Evaluation tests
    score = model.evaluate(xDataTest, yDataTest)
    print('Test score:', score[0])
    print('Test accuracy:', score[1])

    # save model
    model.save('model_M4')
    return 0


# Build the network:
def create_model():
  model = tf.keras.models.Sequential([
    keras.layers.Dense(256, activation='relu', input_shape=(4,)),
    keras.layers.Dropout(0.2),
    keras.layers.Dense(512, activation='relu', input_shape=(64,)),
    keras.layers.Dropout(0.2),
    keras.layers.Dense(256, activation='relu', input_shape=(256,)),
    keras.layers.Dropout(0.2),
    keras.layers.Dense(128, activation='relu', input_shape=(128,)),
    keras.layers.Dropout(0.2),
    keras.layers.Dense(4,)
  ])
  model.summary()
  model.compile(loss=tf.keras.losses.MeanSquaredError(), optimizer='adam', metrics=['accuracy'])

  return model

#Create Callbacks
def create_callbacks():
    cb_list = [keras.callbacks.EarlyStopping(
    monitor='val_loss', min_delta=0, patience=0, verbose=0, mode='auto',
    baseline=None, restore_best_weights=False
)]
    return cb_list

def prepare_data(filename):
    # reading csv file
    dataFrameInput = pd.read_csv(filename) #outputs a dataframe object

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
    #print(DataArray.shape)

    # Strip off header information, i.e. the first 3 cols
    DataArraySlim = DataArray[:, :, 3:]
    #print(DataArraySlim.shape)

    # split in train and test data (ratio 4:1)
    DataTrain = DataArraySlim[:4 * int(DataArraySlim.shape[0] / 5)]
    DataTest = DataArraySlim[4 * int(DataArraySlim.shape[0] / 5):]

    # Split in x (input) and y (output) data
    xDataTrain = DataTrain[:, 0, :]
    yDataTrain = DataTrain[:, 1, :]
    xDataTest = DataTest[:, 0, :]
    yDataTest = DataTest[:, 1, :]

    return (xDataTrain,yDataTrain,xDataTest,yDataTest)


if __name__ == '__main__':
    main()
