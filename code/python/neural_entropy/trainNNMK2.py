### This is a script for the training of the
### Second NN approach

#imports
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import math
from tensorflow.keras.callbacks import EarlyStopping,ModelCheckpoint
import random
import pickle

# Custom Loss
def custom_loss1dMB(u_input, alpha_pred):  # (label,prediciton)
    return 4 * math.pi * tf.math.exp(alpha_pred * np.sqrt(1 / (4 * np.pi))) - alpha_pred * u_input


# Custom Loss
def custom_loss1dMBPrime(): # (label,prediciton)
    def loss(u_input, alpha_pred):
        return 0.5*tf.square(4*math.pi*np.sqrt(1/(4*np.pi))*tf.math.exp(alpha_pred*np.sqrt(1/(4*np.pi))) - u_input)
    return loss

# Build the network:
def create_model():
    # Define the input
    input_ = keras.Input(shape=(1,))

    # Hidden layers
    hidden1 = layers.Dense(4, activation="tanh")(input_)
    hidden2 = layers.Dense(8, activation="tanh")(hidden1)
    hidden3 = layers.Dense(32, activation="tanh")(hidden2)
    hidden4 = layers.Dense(8, activation="tanh")(hidden3)
    hidden5 = layers.Dense(4, activation="tanh")(hidden4)

    # Define the ouput
    output_ = layers.Dense(1)(hidden5)

    # Create the model
    model = keras.Model(inputs=[input_], outputs=[output_])

    model.summary()

    # tf.keras.losses.MeanSquaredError()
    # custom_loss1d
    model.compile(loss=custom_loss1dMBPrime(), optimizer='adam')#, metrics=[custom_loss1dMB, custom_loss1dMBPrime])

    return model

def main():
    
    print("Create Model")
    model = create_model()
	
    print("Create Training Data")
    # build training data and shuffe!
    uTrain = np.arange(0.1, 400, 0.000001)
    random.shuffle(uTrain)
   

    # Create Early Stopping callback
    es = EarlyStopping(monitor='loss', mode='min', min_delta=0.00005, patience=50,
                     verbose=10)  # loss == custom_loss1dMBPrime by model definition
    mc = ModelCheckpoint('saved_model/best_model_1_300.h5', monitor='loss', mode='min', save_best_only=True)

    # Train the model
    print("Train Model")
    history = model.fit(uTrain, uTrain, validation_split=0.3, epochs=1500, batch_size=90000, verbose=1,
                        callbacks=[es, mc])

	
	#save trained model
    print("save model")
    model.save('saved_model/_EntropyLoss_1_300_M_0')
    
    # summarize history for loss
    print("save history")
    with open('saved_model/_EntropyLoss_1_300_M_0_hist.pickle', 'wb') as file_pi:
        pickle.dump(history.history, file_pi)

    # load history
    '''
    history = pickle.load(open('saved_model/_EntropyLoss_1_300_M_0_hist.pickle'), "rb")
    '''
    
    print("Training Sequence successfully finished")
    return 0
   

if __name__ == '__main__':
    main()
