# imports
import tensorflow as tf
import numpy as np


def initialize_network():
    # Load model
    model = tf.keras.models.load_model('saved_model/my_model')

    # Check its architecture
    model.summary()

    return model


# make the network a gobal variable here
model = initialize_network()


def call_network(input):
    inputNP = np.asarray([input])

    predictions = model.predict(inputNP)

    return predictions[0]


def call_networkBatchwise(input):

    #print(input)
    inputNP = np.asarray(input)
    #print(inputNP.shape)
    #print(inputNP)

    predictions = model.predict(inputNP)

    #print(predictions)

    size = predictions.shape[0]*predictions.shape[1]
    test = np.zeros(size)
    for i in  range(0,size):
        test[i] = predictions.flatten(order='C')[i]
    return test

    #return predictions.flatten(order='C')[0]


def main():
    input = [[[0, 1, 2, 3], [2, 3, 4, 5]]]

    # initialize_network()
    print(call_network([0, 1, 2, 3]))
    print("-----")
    print(call_networkBatchwise(input))

    return 0


if __name__ == '__main__':
    main()
