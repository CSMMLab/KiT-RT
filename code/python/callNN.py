#imports
import tensorflow as tf
import numpy as np
import tensorflow.keras.backend as K
import os



def initialize_network():
    # Load model
    model = tf.keras.models.load_model('saved_model/my_model')

    # Check its architecture
    model.summary()

    print("network initizalized")
    return model


# make the network a gobal variable here
model = initialize_network()

def call_network(input):


    #print(tf.__version__)

    inputNP = np.asarray([input])
    #print(os.getcwd())
    #model = tf.keras.models.load_model('saved_model/my_model')
    #model.summary()

   # print("input")
   # print(inputNP)

    predictions = model.predict(inputNP)

    #print(predictions)
    #print(predictions[0])
    #print("python out")

    #print(predictions[0].shape)
    #print(np.zeros(4).shape)

    output = np.zeros(4)
    for i in range(0, 4):
        output[i] = predictions[0, i]

    #print(output)
    return output



def main():
    input = [0,1,2,3]

    #initialize_network()
    print(call_network(input))

    return 0


if __name__ == '__main__':
    main()
