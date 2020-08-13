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
from sklearn.preprocessing import normalize

def main():


    (xDataTrain,yDataTrain,xDataTest,yDataTest) = preprocess_data("trainNN.csv")

    #plot_data(xDataTrain,yDataTrain)

    model = create_model()

    cb_list = create_callbacks()
    # Do the training
    history = model.fit(xDataTrain, yDataTrain, validation_split=0.25, epochs=100, batch_size=500, verbose=1)

    # Evaluation tests
    score = model.evaluate(xDataTest, yDataTest)
    print('Test score:', score[0])
    print('Test accuracy:', score[1])

    # plot training results
    print_output(history)

    # Evaluation tests
    score = model.evaluate(xDataTest, yDataTest)
    print('Test score:', score[0])
    print('Test accuracy:', score[1])

    # save model
    model.save('model_M4')
    return 0


# Build the network:
def create_model():
  ''' mark 1 model
  model = tf.keras.models.Sequential([
    keras.layers.Dense(256, activation='relu', input_shape=(4,)),
    keras.layers.Dropout(0.3),
    keras.layers.Dense(512, activation='relu', input_shape=(64,)),
    keras.layers.Dropout(0.3),
    keras.layers.Dense(256, activation='relu', input_shape=(256,)),
    keras.layers.Dropout(0.3),
    keras.layers.Dense(128, activation='relu', input_shape=(128,)),
    keras.layers.Dropout(0.3),
    keras.layers.Dense(4,)
  ])
  '''
  #leakyRelu = tf.keras.layers.LeakyReLU(alpha=0.1)

  model = tf.keras.models.Sequential([
      keras.layers.Dense(256, activation='sigmoid', input_shape=(4,)),
      keras.layers.Dense(512, activation='sigmoid'),
      keras.layers.Dropout(0.2),
      keras.layers.Dense(256, activation='sigmoid'),
      keras.layers.Dropout(0.2),
      keras.layers.Dense(128, activation='sigmoid'),
      keras.layers.Dense(4, )
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

def preprocess_data(filename):
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

    #Normalize Input
    xDataTrain = normalize(xDataTrain, axis=1, norm='l1')
    xDataTest = normalize(xDataTest, axis=1, norm='l1')

    return (xDataTrain,yDataTrain,xDataTest,yDataTest)

def plot_data(xDataTrain,yDataTrain):
    fig, axs = plt.subplots(4)
    axs[0].plot(xDataTrain[:,0])
    axs[0].set_title('m_0^0')
    axs[0].set_ylabel('value')
    axs[0].set_xlabel('sample id')

    axs[1].plot(xDataTrain[:,1])
    axs[1].set_title('m_1^-1')
    axs[1].set_ylabel('value')
    axs[1].set_xlabel('sample id')

    axs[2].plot(xDataTrain[:,2])
    axs[2].set_title('m_1^0')
    axs[2].set_ylabel('value')
    axs[2].set_xlabel('sample id')

    axs[3].plot(xDataTrain[:,3])
    axs[3].set_title('m_1^1')
    axs[3].set_ylabel('value')
    axs[3].set_xlabel('sample id')
    plt.show()

    fig, axs = plt.subplots(4)
    axs[0].plot(yDataTrain[:, 0])
    axs[0].set_title('m_0^0')
    axs[0].set_ylabel('value')
    axs[0].set_xlabel('sample id')

    axs[1].plot(yDataTrain[:, 1])
    axs[1].set_title('m_1^-1')
    axs[1].set_ylabel('value')
    axs[1].set_xlabel('sample id')

    axs[2].plot(yDataTrain[:, 2])
    axs[2].set_title('m_1^0')
    axs[2].set_ylabel('value')
    axs[2].set_xlabel('sample id')

    axs[3].plot(yDataTrain[:, 3])
    axs[3].set_title('m_1^1')
    axs[3].set_ylabel('value')
    axs[3].set_xlabel('sample id')
    plt.show()

    return 0

def print_output(history):
    fig, axs = plt.subplots(2)
    print(history.history.keys())
    axs[0].plot(history.history['accuracy'])
    axs[0].plot(history.history['val_accuracy'])
    axs[0].set_title('model accuracy')
    axs[0].set_ylabel('accuracy')
    axs[0].set_xlabel('epoch')
    axs[0].legend(['train_acc', 'val_acc'], loc='upper left')
    # axs[0].show()
    # summarize history for loss
    axs[1].plot(history.history['loss'])
    axs[1].plot(history.history['val_loss'])
    axs[1].set_title('model loss')
    axs[1].set_ylabel('loss')
    axs[1].set_xlabel('epoch')
    axs[1].legend(['train_loss', 'val_loss'], loc='upper left')
    # axs[1].show()
    plt.show()
    return 0

if __name__ == '__main__':
    main()
