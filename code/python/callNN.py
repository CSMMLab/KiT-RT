#imports
import tensorflow as tf
import numpy as np
import tensorflow.keras.backend as K
# make the network a gobal variable here

#model = tf.python.keras.engine.sequential.Sequential()
    #tf.keras.Model()

def call_network(input):
    inputNP = np.asarray([input])
    model = tf.keras.models.load_model('saved_model/my_model')

    print("input")
    print(inputNP)

    predictions = model.predict(inputNP)

    print("network called")

    print(predictions)
    print(predictions[0])
    print("python out")

    print(predictions[0].shape)
    print(np.zeros(4).shape)

    output = np.zeros(4)
    for i in range(0, 3):
        output[i] = predictions[0, i]

    return output


def initialize_network():
    # Load model
    model = tf.keras.models.load_model('saved_model/my_model')

    # Check its architecture
    model.summary()

    print("network initizalized")
    return 0


def main():
    input = [0,1,2,3]

    initialize_network()
    call_network(input)

    return 0


if __name__ == '__main__':
    main()
